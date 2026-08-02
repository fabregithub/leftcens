#!/usr/bin/env bash
# =============================================================================
# Launch the D1 definitive validation run(s), DETACHED, so you can close the
# terminal / walk away. Runs each model SEQUENTIALLY (each run already uses
# ~22 of your 24 cores, so running them at once would oversubscribe). Progress
# + ETA stream to per-model logs; a master log tracks the overall sequence.
#
#   ./validation/run_d1.sh                    # tobit then ridge, full grid, 500 reps, M=50
#   MODELS=tobit ./validation/run_d1.sh       # tobit only
#   N_REP=100 M=25 ./validation/run_d1.sh     # a smaller confirmation run
#
# Watch it:   tail -f validation/logs/d1_master_<stamp>.log
#             tail -f validation/logs/d1_tobit_<stamp>.log
# Stop it:    kill $(cat validation/logs/d1_master_<stamp>.pid)   # stops the sequence
# When done:  results in validation/results/d1_<model>_<stamp>.rds (+ _latest per model)
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."                     # repo root

# One BLAS thread per worker: the OpenBLAS pthreads build would otherwise let
# every fork spawn its own thread pool and thrash the machine.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1

# Models to run, in order. Space-separated. Default: the capstone (tobit) then
# the ridge baseline on the same grid.
export MODELS="${MODELS:-tobit ridge}"

export CONFIG="${CONFIG:-full}"
export N_REP="${N_REP:-500}"
export M="${M:-50}"
export ITERS_ALL="${ITERS_ALL:-30}"
export WORKERS="${WORKERS:-22}"
export N_CHUNKS="${N_CHUNKS:-100}"
export BASE_SEED="${BASE_SEED:-20240101}"

mkdir -p validation/results validation/logs
export STAMP="$(date +%Y%m%d-%H%M%S)"
MASTER_LOG="validation/logs/d1_master_${STAMP}.log"
PIDF="validation/logs/d1_master_${STAMP}.pid"

echo "Launching D1 sequence: models=[${MODELS}] config=${CONFIG} N_REP=${N_REP} M=${M} workers=${WORKERS}"
echo "Master log: ${MASTER_LOG}"

# The whole sequence runs inside one detached process so a single kill stops it,
# and it survives closing the terminal.
nohup bash -c '
  echo "[$(date +%H:%M:%S)] D1 sequence starting: models=[${MODELS}]"
  for m in ${MODELS}; do
    export IMP_MODEL="${m}"
    LOG="validation/logs/d1_${m}_${STAMP}.log"
    echo "[$(date +%H:%M:%S)] >>> starting model=${m}  (log: ${LOG})"
    if Rscript validation/run_d1_parallel.R > "${LOG}" 2>&1; then
      echo "[$(date +%H:%M:%S)] <<< finished model=${m}  OK"
    else
      rc=$?
      echo "[$(date +%H:%M:%S)] <<< model=${m} FAILED (exit ${rc}); continuing to next"
    fi
  done
  echo "[$(date +%H:%M:%S)] D1 sequence complete."
' > "${MASTER_LOG}" 2>&1 &

echo $! > "${PIDF}"
echo "PID $(cat "${PIDF}") (saved to ${PIDF})"
echo "Watch:  tail -f ${MASTER_LOG}"
echo "        tail -f validation/logs/d1_tobit_${STAMP}.log"
echo "Stop:   kill \$(cat ${PIDF})"
