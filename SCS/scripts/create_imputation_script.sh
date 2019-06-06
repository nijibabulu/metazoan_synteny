#! /bin/bash

# This script generates a list of commands to make the necessary imputation
# and correlation computations for data imputation. This includes the settings
# from the paper.

for sp in aq ml ta nv sm; do 
    echo Rscript scripts/impute.R --species=$sp --expression=ess --method=none 0 unimputed.$sp 
    echo Rscript scripts/impute.R --species=$sp --expression=ess.norm --method=none 0 unimputed.norm.$sp 

    for t in 5 6 9 17; do
        echo Rscript scripts/impute.R --species=$sp --method=magic $t magic.$sp.$t
    done 

    # the following were not feasible for schmidtea.
    if [ $sp == "sm" ]; then
        continue
    fi
    for k in 5 10 15; do
        echo Rscript scripts/impute.R --species=$sp --method=drimpute $k drimpute.$sp.$k
        echo Rscript scripts/impute.R --species=$sp --method=scimpute $k scimpute.$sp.$k
    done 
done
