#!/bin/bash
# Run DISORT-option experiments listed on stdin as "CASE EXPERIMENT [REPEAT]" lines,
# N at a time, in the controlled environment of validation/tmp/disort_options/env.sh.
#   validation/v24/disort_options/launch.sh N < jobs.txt
export PATH=/usr/bin:/bin
cd /home/g/grainger/project-oraclut
T=validation/tmp/disort_options
mkdir -p $T/logs
xargs -P "${1:-8}" -L 1 bash -c 'env -i bash -c "source '$T'/env.sh; cd /home/g/grainger/project-oraclut; nice -n 19 \$PY validation/v24/disort_options/run_variant.py $0 $1 $2" > '$T'/logs/$0_$1${2:+_r$2}.log 2>&1; echo "$0 $1 $2 exit $?" >> '$T'/logs/status'
