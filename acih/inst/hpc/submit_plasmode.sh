#!/bin/bash
## plasmode study: prepare (waits) -> array -> summary.  Needs app/results from the data analysis.
## Run from the run directory:  bash hpc/submit_plasmode.sh
set -e
export R_LIBS=${ACIH_LIB:-$HOME/Rlib}${R_LIBS:+:$R_LIBS}
mkdir -p logs
sed "s/STEP/prepare/g" hpc/09_plasmode_other.lsf | bsub -K
grep -h "lambda0" $(ls -t logs/pls_prepare.*.err | head -1)
NTASK=$(Rscript acih.R plasmode_ntask)
echo "number of plasmode tasks: $NTASK"
J1=$(sed "s/NTASK/$NTASK/" hpc/08_plasmode_array.lsf | bsub | sed -E 's/.*<([0-9]+)>.*/\1/')
J2=$(sed "s/STEP/summarize/g" hpc/09_plasmode_other.lsf | bsub -w "ended($J1)" | sed -E 's/.*<([0-9]+)>.*/\1/')
echo "submitted: array $J1, summary $J2"
