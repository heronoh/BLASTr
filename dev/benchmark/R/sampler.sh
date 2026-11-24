#!/usr/bin/env bash
# Background sampler for one benchmark run.
# Usage: sampler.sh <run_id> <out_tsv> <interval_s>
# Every <interval_s> seconds writes one line with system and run-level metrics.
# Run processes are found by the env var BLASTR_BENCH_RUN=<run_id>.
# Stops when the file <out_tsv>.stop exists.
# watched_blast_*: blastn of BLASTR_BENCH_WATCH_USER (0 when it is unset).

run_id="$1"
watch_user="${BLASTR_BENCH_WATCH_USER:-}"
out="$2"
interval="${3:-2}"

printf 't\tcpu_busy_jiffies\tcpu_total_jiffies\tsda_read_sectors\tsda_io_ms\tload1\tmem_avail_kb\tcached_kb\trun_nproc\trun_pss_kb\trun_read_bytes\twatched_blast_running\twatched_blast_stopped\trun_cpu_ticks\n' >"$out"

run_pids() {
  local p
  for p in $(pgrep -u "$(id -u)" -f 'blastn|R --no-echo|exec/R' 2>/dev/null); do
    if tr '\0' '\n' <"/proc/$p/environ" 2>/dev/null | grep -qx "BLASTR_BENCH_RUN=$run_id"; then
      echo "$p"
    fi
  done
}

while [[ ! -e "$out.stop" ]]; do
  t=$(date +%s.%N)
  read -r _ u n s idle iow irq sirq st _ </proc/stat
  busy=$((u + n + s + irq + sirq + st))
  total=$((busy + idle + iow))
  read -r rs ioms < <(awk '$3=="sda"{print $6, $13}' /proc/diskstats)
  load1=$(cut -d' ' -f1 /proc/loadavg)
  mavail=$(awk '/^MemAvailable/{print $2}' /proc/meminfo)
  cached=$(awk '/^Cached/{print $2}' /proc/meminfo)
  pids=$(run_pids)
  np=0; pss=0; rb=0; ct=0
  for p in $pids; do
    np=$((np + 1))
    v=$(awk '/^Pss:/{print $2}' "/proc/$p/smaps_rollup" 2>/dev/null); pss=$((pss + ${v:-0}))
    v=$(awk '/^read_bytes/{print $2}' "/proc/$p/io" 2>/dev/null); rb=$((rb + ${v:-0}))
    # own utime+stime only (fields 14-15); per-pid lines go to <out>.pids so
    # the analysis can take the max per pid and sum (no double counting)
    v=$(sed 's/.*) //' "/proc/$p/stat" 2>/dev/null | awk '{print $12+$13}'); ct=$((ct + ${v:-0}))
    printf '%s\t%s\t%s\n' "$p" "$(cat /proc/$p/comm 2>/dev/null)" "${v:-0}" >>"$out.pids"
  done
  hr=0; hs=0
  if [[ -n "$watch_user" ]]; then
    hr=$(ps -u "$watch_user" -o stat=,comm= 2>/dev/null | awk '$2=="blastn" && $1 !~ /^T/' | wc -l)
    hs=$(ps -u "$watch_user" -o stat=,comm= 2>/dev/null | awk '$2=="blastn" && $1 ~ /^T/' | wc -l)
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$t" "$busy" "$total" "$rs" "$ioms" "$load1" "$mavail" "$cached" "$np" "$pss" "$rb" "$hr" "$hs" "$ct" >>"$out"
  sleep "$interval"
done
