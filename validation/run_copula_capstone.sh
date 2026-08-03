#!/usr/bin/env bash
# =============================================================================
# Launch the copula capstone DETACHED (walk away / close the terminal).
# Definitive high-replication validation of imp_model="copula" on the full grid,
# using the recommended hybrid MI (plug-in point + drawn-margin variance).
# Same grid/seeds as the tobit D1, so results pair directly with
# validation/results/d1_tobit_latest.rds.
#
#   ./validation/run_copula_capstone.sh                 # N_REP=300 M=30 (~12 h)
#   N_REP=500 M=50 ./validation/run_copula_capstone.sh  # exact D1 parity (~30 h)
#
# Watch:  tail -f validation/logs/copula_capstone_<stamp>.log
# Stop:   kill $(cat validation/logs/copula_capstone_<stamp>.pid)
# Done:   validation/results/copula_capstone_latest.rds (+ timestamped, csv, summary)
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export CONFIG="${CONFIG:-full}"
export N_REP="${N_REP:-300}"
export M="${M:-30}"
export ITERS_ALL="${ITERS_ALL:-30}"
export WORKERS="${WORKERS:-22}"
export N_CHUNKS="${N_CHUNKS:-100}"
export BASE_SEED="${BASE_SEED:-20240101}"

mkdir -p validation/results validation/logs
STAMP="$(date +%Y%m%d-%H%M%S)"
LOG="validation/logs/copula_capstone_${STAMP}.log"
PIDF="validation/logs/copula_capstone_${STAMP}.pid"

echo "Launching copula capstone: config=${CONFIG} N_REP=${N_REP} M=${M} workers=${WORKERS}"
nohup Rscript validation/run_copula_capstone.R > "${LOG}" 2>&1 &
echo $! > "${PIDF}"
echo "PID $(cat "${PIDF}")  |  log: ${LOG}"
echo "Watch: tail -f ${LOG}"
echo "Stop:  kill \$(cat ${PIDF})"
