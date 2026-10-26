#!/bin/bash
## prepare -> simulation array -> summary, chained with LSF dependencies.
## Run from the run directory:  bash hpc/submit_all.sh
## If out/designs.rds already exists (prepare done), only the array and the summary are submitted.
set -e
mkdir -p logs out results
export R_LIBS=${ACIH_LIB:-$HOME/Rlib}${R_LIBS:+:$R_LIBS}
if [ ! -f out/designs.rds ]; then
  echo "designs not prepared yet: submitting prepare and waiting for it (bsub -K)..."
  bsub -K < hpc/00_prepare.lsf
fi
NTASK=$(Rscript acih.R ntask)
echo "number of simulation tasks: $NTASK"
J1=$(sed "s/NTASK/$NTASK/" hpc/01_sim_array.lsf | bsub | sed -E 's/.*<([0-9]+)>.*/\1/')
J2=$(bsub -w "ended($J1)" < hpc/02_summarize.lsf | sed -E 's/.*<([0-9]+)>.*/\1/')
echo "submitted: array $J1, summary $J2"
