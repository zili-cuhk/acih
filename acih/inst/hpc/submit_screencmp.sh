#!/bin/bash
## screening-criterion check: array -> summary.  Needs out/designs.rds from the simulation set-up.
## Run from the run directory:  bash hpc/submit_screencmp.sh
set -e
export R_LIBS=${ACIH_LIB:-$HOME/Rlib}${R_LIBS:+:$R_LIBS}
mkdir -p logs
NTASK=$(Rscript acih.R screencmp_ntask)
echo "number of screening-check tasks: $NTASK"
J1=$(sed "s/NTASK/$NTASK/" hpc/10_screencmp_array.lsf | bsub | sed -E 's/.*<([0-9]+)>.*/\1/')
J2=$(bsub -w "ended($J1)" < hpc/11_screencmp_summarize.lsf | sed -E 's/.*<([0-9]+)>.*/\1/')
echo "submitted: array $J1, summary $J2"
