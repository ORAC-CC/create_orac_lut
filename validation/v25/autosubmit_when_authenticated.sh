#!/usr/bin/env bash
# Submit the definitive V25 production set as soon as a valid Kerberos ticket is
# available (validation helper).  The established submission wrapper hops to
# atmlxint7 over ssh, which needs the user's ticket (sshd must read
# ~/.ssh/authorized_keys on the Kerberos-protected home).  This script polls
# klist every two minutes for up to seven days, then runs
#    scripts/submit_v25_cloud_luts.sh all
# once, verifies the submission and exits.  Log: validation/tmp/v25/logs/autosubmit.log
#   nohup validation/v25/autosubmit_when_authenticated.sh &      (started 2026-10-06)
#   pkill -f autosubmit_when_authenticated                          (to cancel)
set -u
export PATH=/usr/bin:/bin
cd /home/g/grainger/project-oraclut
LOG=validation/tmp/v25/logs/autosubmit.log
EXPECTED_REV=990d970e19e1411db44121153b110c8ff5cfae09
echo "$(date -u +%FT%TZ) waiting for a valid Kerberos ticket (cache ${KRB5CCNAME:-default})" >> "$LOG"
for ((i = 0; i < 5040; i++)); do
    if klist -s 2>/dev/null; then
        echo "$(date -u +%FT%TZ) ticket valid; checking HEAD" >> "$LOG"
        [[ "$(git rev-parse HEAD)" == "$EXPECTED_REV" ]] || { echo "HEAD is not $EXPECTED_REV; not submitting" >> "$LOG"; exit 2; }
        already=$(awk -F'\t' -v rev="$EXPECTED_REV" 'NR>1 && $3==rev' validation/v25/production_submissions.tsv | wc -l)
        if [[ "$already" -gt 0 ]]; then echo "$(date -u +%FT%TZ) manifest already has $already rows at $EXPECTED_REV (submitted elsewhere); not submitting" >> "$LOG"; exit 0; fi
        if ! ssh -o BatchMode=yes -o ConnectTimeout=20 atmlxint7 true >> "$LOG" 2>&1; then
            echo "$(date -u +%FT%TZ) ssh to atmlxint7 still failing; retrying" >> "$LOG"; sleep 120; continue
        fi
        echo "$(date -u +%FT%TZ) submitting" >> "$LOG"
        SBATCH_PARTITION=shared,priority-eodg bash scripts/submit_v25_cloud_luts.sh all >> "$LOG" 2>&1
        status=$?
        echo "$(date -u +%FT%TZ) submit script exit $status" >> "$LOG"
        n=$(awk -F'\t' -v rev="$EXPECTED_REV" 'NR>1 && $3==rev' validation/v25/production_submissions.tsv | wc -l)
        echo "$(date -u +%FT%TZ) manifest rows at $EXPECTED_REV: $n" >> "$LOG"
        bash validation/v25/monitor_jobs.sh >> "$LOG" 2>&1
        exit $status
    fi
    sleep 120
done
echo "$(date -u +%FT%TZ) gave up after seven days" >> "$LOG"
exit 3
