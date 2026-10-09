#!/usr/bin/env bash
# Show the state of the V25 production jobs recorded in the submission manifest
# (validation only).  Queries SLURM on atmlxint7 over ssh (the caller's ssh agent
# or Kerberos ticket) and lists, per job, the queue state or the job log's
# completion line and the product's presence.
#   validation/v25/monitor_jobs.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
MANIFEST=validation/v25/production_submissions.tsv
[[ -s "$MANIFEST" ]] || { echo "no submissions recorded in $MANIFEST"; exit 1; }
ids="$(tail -n +2 "$MANIFEST" | cut -f4 | paste -sd, -)"
queue="$(ssh -o BatchMode=yes atmlxint7 "squeue -h -j $ids -o '%i %T %M %L %N' 2>/dev/null" || true)"
printf '%-8s %-26s %-10s %-24s %s\n' JOB NAME STATE TIME/LEFT PRODUCT
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r when host rev id name run product; do
    state="$(awk -v j="$id" '$1 == j {print $2" "$3"/"$4}' <<< "$queue")"
    if [[ -z "$state" ]]; then
        log="$(ls validation/slurm/${name}_${id}_*.out 2>/dev/null | head -1)"
        if [[ -n "$log" ]] && grep -q "ORAC LUT SLURM JOB COMPLETED" "$log"; then state="COMPLETED"
        elif [[ -n "$log" ]]; then state="ENDED (no completion line: see $log)"
        else state="NOT IN QUEUE (no log)"; fi
    fi
    present="missing"; [[ -s "$product" ]] && present="present ($(du -h "$product" | cut -f1))"
    printf '%-8s %-26s %-10s %-24s %s\n' "$id" "$name" "${state%% *}" "${state#* }" "$present"
done
