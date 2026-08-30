#!/bin/bash
#SBATCH --nodes=1

set -euo pipefail

campaign_root=${OQP_BA_CAMPAIGN_ROOT:?set OQP_BA_CAMPAIGN_ROOT}
threads=${OQP_OMP_THREADS:-8}
indices_csv=${OQP_PACK_INDICES:?set comma-separated trajectory indices}
worker="$campaign_root/scripts/run_array.sbatch"
IFS=, read -ra indices <<< "$indices_csv"

declare -A allowed_cpu=()
allowed_list=$(awk '/Cpus_allowed_list/ {print $2}' /proc/self/status)
IFS=, read -ra allowed_parts <<< "$allowed_list"
for part in "${allowed_parts[@]}"; do
  if [[ "$part" == *-* ]]; then
    for cpu in $(seq "${part%-*}" "${part#*-}"); do allowed_cpu[$cpu]=1; done
  else
    allowed_cpu[$part]=1
  fi
done
declare -A seen_core=()
physical_cpus=()
while IFS=, read -r cpu core; do
  [[ "$cpu" == \#* ]] && continue
  if [[ -n ${allowed_cpu[$cpu]:-} && -z ${seen_core[$core]:-} ]]; then
    physical_cpus+=("$cpu")
    seen_core[$core]=1
  fi
done < <(lscpu -p=CPU,CORE)
test "${#physical_cpus[@]}" -ge $((${#indices[@]} * threads))

pids=()
for slot in "${!indices[@]}"; do
  first=$((slot * threads))
  cpu_list=$(IFS=,; echo "${physical_cpus[*]:first:threads}")
  env SLURM_ARRAY_TASK_ID="${indices[$slot]}" SLURM_CPUS_PER_TASK="$threads" \
    taskset -c "$cpu_list" bash "$worker" &
  pids+=("$!")
done
status=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then status=1; fi
done
exit "$status"
