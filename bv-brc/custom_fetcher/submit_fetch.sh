#!/usr/bin/env bash
# Submit the whole accession -> genome download as a Slurm workflow.
#   ./submit_fetch.sh accessions.txt OUTDIR [NSHARDS=20] [INCLUDE=genome]
# Re-running the same command is safe: finished genomes are skipped, only failures/pending are retried,
# and resolution.tsv is NOT regenerated (so replacement decisions stay frozen) unless FORCE_RESOLVE=1.
#
# Cluster-specific knobs (export before calling; none are assumed - check your site docs):
#   SBATCH_EXTRA   extra sbatch args, e.g. "--partition=<name> --account=<acct>"
#   PY_SETUP       shell snippet making python3 (>=3.8) available, e.g. "module load python"
#   ARRAY_CONC     max simultaneous array tasks (default 4)
#   NCBI_API_KEY   optional; raises NCBI rate limit 5 -> 10 req/s (propagated to jobs)
set -euo pipefail
ACC=$(cd "$(dirname "${1:?accessions.txt}")" && pwd)/$(basename "$1"); mkdir -p "${2:?outdir}"; OUT=$(cd "$2" && pwd)
NSHARDS=${3:-20}; INCLUDE=${4:-genome}
ARRAY_CONC=${ARRAY_CONC:-4}; SBATCH_EXTRA=${SBATCH_EXTRA:-}; PY_SETUP=${PY_SETUP:-:}
SCRIPT=$(cd "$(dirname "$0")" && pwd)/fetch_genomes.py
mkdir -p "$OUT/logs"
# NCBI limit is per IP/key: split the budget across concurrently running shards
TOTAL_RPS=$([ -n "${NCBI_API_KEY:-}" ] && echo 8 || echo 4)
RPS=$(awk -v t="$TOTAL_RPS" -v c="$ARRAY_CONC" 'BEGIN{printf "%.2f", t/c}')
COMMON="--parsable --export=ALL $SBATCH_EXTRA --chdir=$OUT"

DEP=""
if [[ ! -s "$OUT/resolution.tsv" || "${FORCE_RESOLVE:-0}" == "1" ]]; then
  J1=$(sbatch $COMMON --job-name=acc_resolve --time=02:00:00 --mem=2G --cpus-per-task=1 \
       --output="$OUT/logs/resolve_%j.log" --wrap="$PY_SETUP; python3 $SCRIPT resolve $ACC $OUT")
  DEP="--dependency=afterok:$J1"; echo "resolve job: $J1"
fi
J2=$(sbatch $COMMON $DEP --job-name=acc_dl --array=0-$((NSHARDS-1))%$ARRAY_CONC --time=12:00:00 \
     --mem=4G --cpus-per-task=2 --output="$OUT/logs/dl_%A_%a.log" \
     --wrap="$PY_SETUP; python3 $SCRIPT download $OUT --nshards $NSHARDS --shard \$SLURM_ARRAY_TASK_ID --threads 2 --rps $RPS --include $INCLUDE")
echo "download array: $J2 ($NSHARDS shards, $ARRAY_CONC concurrent, $RPS req/s each)"
J3=$(sbatch $COMMON --dependency=afterany:$J2 --job-name=acc_merge --time=00:30:00 --mem=2G \
     --output="$OUT/logs/merge_%j.log" --wrap="$PY_SETUP; python3 $SCRIPT merge $OUT")
echo "merge job: $J3 -> $OUT/manifest.tsv, $OUT/summary.txt (merge exits 1 if anything is failed/pending: just re-run this script)"
