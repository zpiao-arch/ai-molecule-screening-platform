#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULT="${HTS_RESULT:-${SCRIPT_DIR}/hts_closed_loop_results.json}"
LOG="${HTS_LOG:-${SCRIPT_DIR}/hts_loop.log}"
echo "monitor started at $(date)"
for i in $(seq 1 120); do
  if [ -f "$RESULT" ]; then
    echo "RESULT_READY at iter=$i ($(date))"
    break
  fi
  if ! pgrep -f "hts_closed_loop.py" >/dev/null 2>&1; then
    echo "PROCESS_ENDED_NO_RESULT at iter=$i ($(date))"
    break
  fi
  sleep 30
done
echo "monitor exit"
