#!/usr/bin/env bash
# =============================================================================
# Launch the E8 m-stability study DETACHED (walk away / close the terminal).
# How many imputations (m) are needed -- stability, FMI, and the effect of
# sample size n -- for copula vs tobit. Manuscript config by default.
#
#   ./validation/run_m_stability.sh                 # REPS=100, NS to 4800 (~4.5 h)
#   NS=75,300,1200 REPS=20 ./validation/run_m_stability.sh   # quick look
#
# Run this AFTER the copula capstone has finished (it wants the whole machine).
# Watch:  tail -f validation/logs/m_stability_<stamp>.log
# Stop:   kill $(cat validation/logs/m_stability_<stamp>.pid)
# Done:   validation/results/effect_of_m_stability.rds (+ *_summary.csv, *_adaptive.csv)
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export REPS="${REPS:-100}"
export NS="${NS:-75,150,300,600,1200,2400,4800}"
export SKEWS="${SKEWS:-0,0.75}"
export MODELS="${MODELS:-tobit,copula}"
export M_MAX="${M_MAX:-100}"
export ND="${ND:-0.35}"
export ITERS="${ITERS:-20}"
export WORKERS="${WORKERS:-22}"

mkdir -p validation/results validation/logs
STAMP="$(date +%Y%m%d-%H%M%S)"
LOG="validation/logs/m_stability_${STAMP}.log"
PIDF="validation/logs/m_stability_${STAMP}.pid"

echo "Launching E8 m-stability: REPS=${REPS} NS=${NS} MODELS=${MODELS} workers=${WORKERS}"
nohup Rscript validation/effect_of_m_stability.R > "${LOG}" 2>&1 &
echo $! > "${PIDF}"
echo "PID $(cat "${PIDF}")  |  log: ${LOG}"
echo "Watch: tail -f ${LOG}"
echo "Stop:  kill \$(cat ${PIDF})"
