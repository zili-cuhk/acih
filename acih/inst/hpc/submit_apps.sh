#!/bin/bash
## data applications: prepare (waits) -> array over cohorts x 100 splits -> summary.
## Run from the run directory:  bash hpc/submit_apps.sh
set -e
mkdir -p logs app/out app/results
export R_LIBS=${ACIH_LIB:-$HOME/Rlib}${R_LIBS:+:$R_LIBS}
echo "preprocessing (bsub -K waits for it)..."
bsub -K < hpc/05_app_prepare.lsf
tail -n 30 $(ls -t logs/app_prepare.*.err | head -1)
NTASK=$(Rscript acih.R app_ntask)
echo "number of application tasks: $NTASK"
J1=$(sed "s/NTASK/$NTASK/" hpc/06_app_array.lsf | bsub | sed -E 's/.*<([0-9]+)>.*/\1/')
J2=$(bsub -w "ended($J1)" < hpc/07_app_summarize.lsf | sed -E 's/.*<([0-9]+)>.*/\1/')
echo "submitted: array $J1, summary $J2"
