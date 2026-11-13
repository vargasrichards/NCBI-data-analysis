# NCBI assembly fetcher (GCA_/GCF_ list -> genomes), Slurm-ready

Stdlib-only Python (>=3.8) against the NCBI Datasets v2 REST API. No conda/containers needed on the cluster.

## Run
    export NCBI_API_KEY=...            # optional, doubles the rate limit
    export SBATCH_EXTRA="--partition=<yours> --account=<yours>"   # site-specific, not assumed
    export PY_SETUP="module load python"                          # if python3 is not on PATH
    ./submit_fetch.sh accessions.txt results_dir 20 genome,gff3
Submits: resolve (1 job) -> download (array of N shards, ARRAY_CONC at a time) -> merge (runs even if some shards fail).
Re-run the same command to retry only failures. `resolution.tsv` is frozen after the first run (FORCE_RESOLVE=1 to redo).

## Input
One accession per line (first column used; `#` comments, blank lines, lower-case, duplicates are tolerated). Versions are required (GCF_000006785.2).

## Retired / superseded entries (resolution.tsv `status`)
| status | meaning | action |
|---|---|---|
| ok | live record at NCBI | downloaded as requested |
| replaced | requested version gone/suppressed; NCBI revision history gives a live successor (`newer_version_same_db`) or the GCA<->GCF counterpart (`other_database_counterpart`) | downloaded; file named by the REQUESTED accession, true accession in `resolved`, NCBI's `identical` flag in `revision_identical_flag` |
| retired | NCBI record exists as suppressed, no live successor | not downloaded |
| not_found | NCBI has no record/successor (withdrawn or typo) | not downloaded |
| invalid | not GCA_/GCF_ or missing version | not downloaded |
`--replacement none` (resolve step) turns off substitution if you need exactly what BV-BRC analysed. Check `replaced` rows before using them: a newer revision can be a different assembly.

## Outputs
`genomes/<requested>.fna.gz` (+ `.gff.gz` etc. with --include; GFF is absent for un-annotated GenBank assemblies and is left blank, not an error),
`manifest.tsv` (every input row, status, resolved accession, path, sha256), `summary.txt`, `provenance.json` (date, input sha256, script sha256), `logs/`.
Reproducibility: gzip written with mtime=0 (byte-identical for identical content); per-file md5 from NCBI's zip is verified before writing; writes are atomic.

## Tests
`python test_classify.py` (offline logic); `python fetch_genomes.py resolve example_accessions.txt out && python fetch_genomes.py download out && python fetch_genomes.py merge out`.
