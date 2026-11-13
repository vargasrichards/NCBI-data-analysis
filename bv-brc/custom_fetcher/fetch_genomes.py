#!/usr/bin/env python3
"""
fetch_genomes.py - reproducible, resumable download of NCBI assemblies (GCA_/GCF_)
from a plain-text accession list, with handling of retired / superseded entries.

Stdlib only (Python >= 3.8). Talks to the NCBI Datasets v2 REST API.

Sub-commands (run in this order; each is idempotent and safe to re-run):
  resolve   accessions.txt -> OUT/resolution.tsv   (current / replaced / suppressed / not_found / invalid)
  download  OUT/resolution.tsv -> OUT/genomes/<requested_acc>.fna.gz [+ .gff.gz ...]  (shardable for Slurm arrays)
  merge     per-shard status files -> OUT/manifest.tsv + OUT/summary.txt

Optional: export NCBI_API_KEY=... (raises rate limit from 5 to 10 requests/s).
"""
import argparse, csv, datetime as dt, gzip, hashlib, http.client, io, json, os, random, re, shutil
import sys, threading, time, urllib.error, urllib.request, zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed

API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
ACC_RE = re.compile(r"^(GC[AF])_(\d{9})(?:\.(\d+))?$")
INCLUDE_MAP = {  # CLI-style name -> (API annotation type, filename suffix in output)
    "genome": ("GENOME_FASTA", ".fna.gz"),
    "gff3": ("GENOME_GFF", ".gff.gz"),
    "gbff": ("GENOME_GBFF", ".gbff.gz"),
    "gtf": ("GENOME_GTF", ".gtf.gz"),
    "cds": ("CDS_FASTA", ".cds.fna.gz"),
    "rna": ("RNA_FASTA", ".rna.fna.gz"),
    "protein": ("PROT_FASTA", ".faa.gz"),
}
RES_COLS = ["requested", "status", "resolved", "resolution", "note", "paired_accession",
            "assembly_name", "organism", "assembly_level", "release_date", "revision_identical_flag"]

# ---------------------------------------------------------------- HTTP layer
class RateLimiter:
    """Thread-safe minimum spacing between requests (per process)."""
    def __init__(self, rps):
        self.gap, self.lock, self.next = 1.0 / rps, threading.Lock(), 0.0
    def wait(self):
        with self.lock:
            now = time.monotonic()
            t = max(now, self.next)
            self.next = t + self.gap
        if t > now:
            time.sleep(t - now)

LIMITER = RateLimiter(3.0)
UA = "ncbi-fetch-genomes/1.0 (reproducible assembly download script)"

def http_get(url, retries=6, timeout=180, binary=False):
    """GET with exponential backoff + jitter; honours Retry-After; returns (status, body).
    4xx other than 429 are returned (not retried); persistent 5xx/network errors raise."""
    headers = {"User-Agent": UA, "Accept": "application/json"}
    if os.environ.get("NCBI_API_KEY"):
        headers["api-key"] = os.environ["NCBI_API_KEY"]
    last = None
    for attempt in range(retries):
        LIMITER.wait()
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                body = r.read()
                return r.status, (body if binary else json.loads(body or b"{}"))
        except urllib.error.HTTPError as e:
            body = e.read()
            if e.code == 429 or e.code >= 500:
                ra = e.headers.get("Retry-After")
                delay = float(ra) if ra and ra.isdigit() else min(2 ** attempt * 2, 120)
                last = f"HTTP {e.code}"
            else:
                try:
                    return e.code, (body if binary else json.loads(body or b"{}"))
                except Exception:
                    return e.code, body
        except (urllib.error.URLError, http.client.HTTPException, TimeoutError, ConnectionError, OSError, json.JSONDecodeError) as e:
            delay, last = min(2 ** attempt * 2, 120), repr(e)
        time.sleep(delay + random.uniform(0, 1.5))
    raise RuntimeError(f"giving up on {url[:120]} after {retries} attempts: {last}")

# ---------------------------------------------------------------- helpers
def read_accessions(path):
    """Normalise (strip, upper-case, drop comments/blank, take first column), de-duplicate keeping order."""
    seen, out = set(), []
    for line in open(path, encoding="utf-8-sig"):
        line = line.split("#")[0].strip()
        if not line:
            continue
        a = re.split(r"[\s,;\t]+", line)[0].strip().upper()
        if a not in seen:
            seen.add(a); out.append(a)
    return out

def chunks(xs, n):
    for i in range(0, len(xs), n):
        yield xs[i:i + n]

def fetch_reports(accs):
    """accession -> dataset_report dict for those NCBI knows as live records. Missing accession => absent."""
    out = {}
    batches = list(chunks(sorted(set(accs)), 50))
    for n, batch in enumerate(batches, 1):
        if n % 20 == 0:
            print(f"[resolve] dataset_report batch {n}/{len(batches)}", file=sys.stderr, flush=True)
        st, j = http_get(f"{API}/genome/accession/{','.join(batch)}/dataset_report")
        if st != 200:
            raise RuntimeError(f"dataset_report HTTP {st}: {str(j)[:200]}")
        for r in j.get("reports", []):
            out[r["accession"]] = r
    return out

def revisions(acc):
    st, j = http_get(f"{API}/genome/accession/{acc}/revision_history")
    return j.get("assembly_revisions", []) if st == 200 and isinstance(j, dict) else []

def meta(rep):
    ai = rep.get("assembly_info", {})
    return dict(paired_accession=rep.get("paired_accession", ""),
                assembly_name=ai.get("assembly_name", ""),
                organism=rep.get("organism", {}).get("organism_name", ""),
                assembly_level=ai.get("assembly_level", ""),
                release_date=ai.get("release_date", ""))

def is_live(rep):
    return rep is not None and rep.get("assembly_info", {}).get("assembly_status", "current") == "current"

# ---------------------------------------------------------------- resolve
def classify(requested, reports, rev_lookup, allow_replace=True):
    """Pure function (unit-testable). Returns a row dict for resolution.tsv."""
    row = {c: "" for c in RES_COLS}
    row["requested"] = requested
    m = ACC_RE.match(requested)
    if not m:
        row.update(status="invalid", note="not a GCA_/GCF_ accession")
        return row
    if not m.group(3):
        row.update(status="invalid", note="no version suffix (e.g. .1); BV-BRC accessions are normally versioned")
        return row
    rep = reports.get(requested)
    if rep is not None:
        row.update(meta(rep))
        astat = rep.get("assembly_info", {}).get("assembly_status", "current")
        if astat == "current":
            row.update(status="ok", resolved=requested, resolution="exact")
            return row
        row["note"] = f"NCBI assembly_status={astat}; "  # suppressed / replaced etc: look for a live successor below
    if not allow_replace:
        row.update(status="retired", note=row["note"] + "not live at NCBI; replacement lookup disabled")
        return row
    # retired / superseded / suppressed: walk revision history, newest first, prefer same database
    prefix = m.group(1)
    same, other = [], []
    for r in rev_lookup(requested):
        flag = r.get("identical")
        for key, bucket in (("genbank_accession", "GCA"), ("refseq_accession", "GCF")):
            cand = r.get(key)
            if cand and cand != requested:
                (same if bucket == prefix else other).append((cand, flag))
    for cand, flag in same + other:
        crep = reports.get(cand)
        if is_live(crep):
            kind = "newer_version_same_db" if cand[:3] == prefix else "other_database_counterpart"
            row.update(meta(crep))
            row.update(status="replaced", resolved=cand, resolution=kind, revision_identical_flag=str(flag),
                       note=row["note"] + f"{requested} retired; using {cand}")
            return row
    row.update(status="not_found" if not row["note"] else "retired",
               note=row["note"] + "no live record or live successor at NCBI (retired/withdrawn/typo)")
    return row

def cmd_resolve(a):
    accs = read_accessions(a.accessions)
    os.makedirs(a.outdir, exist_ok=True)
    print(f"[resolve] {len(accs)} unique accessions", file=sys.stderr)
    valid = [x for x in accs if ACC_RE.match(x) and ACC_RE.match(x).group(3)]
    reports = fetch_reports(valid)
    # second pass: candidate replacements for anything not live
    cache = {}
    def rev_lookup(acc):
        if acc not in cache:
            cache[acc] = revisions(acc)
            if len(cache) % 50 == 0:
                print(f"[resolve] revision lookups done: {len(cache)}", file=sys.stderr, flush=True)
        return cache[acc]
    need = [x for x in valid if not is_live(reports.get(x))]
    print(f"[resolve] {len(valid) - len(need)} live, {len(need)} need replacement lookup", file=sys.stderr, flush=True)
    cands = set()
    if a.replacement != "none":
        for x in need:
            for r in rev_lookup(x):
                cands.update(v for v in (r.get("genbank_accession"), r.get("refseq_accession")) if v)
        cands -= set(reports)
        reports.update(fetch_reports(sorted(cands)))
    rows = [classify(x, reports, rev_lookup, allow_replace=(a.replacement != "none")) for x in accs]
    tmp = os.path.join(a.outdir, "resolution.tsv.tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, RES_COLS, delimiter="\t"); w.writeheader(); w.writerows(rows)
    os.replace(tmp, os.path.join(a.outdir, "resolution.tsv"))
    prov = dict(date_utc=dt.datetime.utcnow().isoformat() + "Z", api=API, input=os.path.abspath(a.accessions),
                input_sha256=hashlib.sha256(open(a.accessions, "rb").read()).hexdigest(),
                replacement_policy=a.replacement, script_sha256=hashlib.sha256(open(__file__, "rb").read()).hexdigest())
    json.dump(prov, open(os.path.join(a.outdir, "provenance.json"), "w"), indent=2)
    from collections import Counter
    c = Counter((r["status"], r["resolution"]) for r in rows)
    for k, v in sorted(c.items()):
        print(f"[resolve] {k[0]:10s} {k[1]:28s} {v}", file=sys.stderr)

# ---------------------------------------------------------------- download
def md5_of(fp, bufsize=1 << 20):
    h = hashlib.md5()
    for b in iter(lambda: fp.read(bufsize), b""):
        h.update(b)
    return h.hexdigest()

def download_one(resolved, include, outdir):
    """Fetch zip for one accession, verify, write gz outputs atomically. Returns dict of output paths + sha256."""
    types = "&".join(f"include_annotation_type={INCLUDE_MAP[i][0]}" for i in include)
    st, data = http_get(f"{API}/genome/accession/{resolved}/download?{types}&hydrated=FULLY_HYDRATED",
                        binary=True, timeout=600)
    if st != 200:
        raise RuntimeError(f"download HTTP {st}: {data[:200]!r}")
    z = zipfile.ZipFile(io.BytesIO(data))
    bad = z.testzip()
    if bad:
        raise RuntimeError(f"corrupt zip member {bad}")
    md5s = {}
    for line in z.read("md5sum.txt").decode().splitlines():
        h, _, p = line.partition("  "); md5s[p.strip()] = h
    outs = {}
    for inc in include:
        suffix = INCLUDE_MAP[inc][1]
        pats = {"genome": lambda n: n.endswith("_genomic.fna"), "gff3": lambda n: n.endswith("genomic.gff"),
                "gbff": lambda n: n.endswith("genomic.gbff"), "gtf": lambda n: n.endswith("genomic.gtf"),
                "cds": lambda n: n.endswith("cds_from_genomic.fna"), "rna": lambda n: n.endswith("rna.fna"),
                "protein": lambda n: n.endswith("protein.faa")}[inc]
        members = [n for n in z.namelist() if f"/{resolved}/" in n and pats(n)]
        if not members:
            if inc == "genome":
                raise RuntimeError("no genomic FASTA in zip")
            outs[inc] = None            # e.g. GFF absent for GenBank-only assemblies: not an error
            continue
        member = members[0]
        with z.open(member) as fh:
            raw = fh.read()
        if hashlib.md5(raw).hexdigest() != md5s.get(member):
            raise RuntimeError(f"md5 mismatch for {member}")
        if inc == "genome" and not raw.startswith(b">"):
            raise RuntimeError("FASTA does not start with '>'")
        outs[inc] = (member, raw, suffix)
    return outs

def write_outputs(requested, outs, outdir):
    res = {}
    for inc, v in outs.items():
        if v is None:
            res[inc] = ""; continue
        member, raw, suffix = v
        final = os.path.join(outdir, "genomes", requested + suffix)
        tmp = final + f".tmp{os.getpid()}"
        with open(tmp, "wb") as fh:
            with gzip.GzipFile(fileobj=fh, mode="wb", mtime=0, compresslevel=6) as g:   # mtime=0 => deterministic bytes
                g.write(raw)
        os.replace(tmp, final)
        res[inc] = hashlib.sha256(open(final, "rb").read()).hexdigest()
    return res

def cmd_download(a):
    global LIMITER
    LIMITER = RateLimiter(a.rps)
    include = [i.strip() for i in a.include.split(",") if i.strip()]
    for i in include:
        if i not in INCLUDE_MAP: sys.exit(f"unknown --include {i}; choose from {list(INCLUDE_MAP)}")
    rows = list(csv.DictReader(open(os.path.join(a.outdir, "resolution.tsv")), delimiter="\t"))
    todo = [r for r in rows if r["status"] in ("ok", "replaced")]
    todo = [r for i, r in enumerate(todo) if i % a.nshards == a.shard]            # deterministic sharding
    os.makedirs(os.path.join(a.outdir, "genomes"), exist_ok=True)
    os.makedirs(os.path.join(a.outdir, "status"), exist_ok=True)
    sfile = os.path.join(a.outdir, "status", f"shard_{a.shard:04d}.tsv")
    done = {}
    if os.path.exists(sfile):                                                      # resume: keep previous successes
        for r in csv.DictReader(open(sfile), delimiter="\t"):
            if r["download"] == "ok" and os.path.exists(os.path.join(a.outdir, "genomes", r["requested"] + ".fna.gz")):
                done[r["requested"]] = r
    cols = ["requested", "resolved", "download", "error", "genome_sha256"] + [f"{i}_sha256" for i in include if i != "genome"]
    print(f"[download] shard {a.shard}/{a.nshards}: {len(todo)} assemblies, {len(done)} already done", file=sys.stderr)
    lock = threading.Lock()
    with open(sfile + ".new", "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore"); w.writeheader()
        for r in done.values(): w.writerow(r)
        fh.flush()
        def work(r):
            rec = {"requested": r["requested"], "resolved": r["resolved"]}
            try:
                outs = download_one(r["resolved"], include, a.outdir)
                sh = write_outputs(r["requested"], outs, a.outdir)
                rec.update(download="ok", error="", genome_sha256=sh.get("genome", ""))
                for i in include:
                    if i != "genome": rec[f"{i}_sha256"] = sh.get(i, "")
            except Exception as e:
                rec.update(download="failed", error=str(e)[:300])
            return rec
        pending = [r for r in todo if r["requested"] not in done]
        with ThreadPoolExecutor(a.threads) as ex:
            futs = [ex.submit(work, r) for r in pending]
            for n, f in enumerate(as_completed(futs), 1):
                rec = f.result()
                with lock:
                    w.writerow(rec); fh.flush()
                if n % 50 == 0 or rec["download"] != "ok":
                    print(f"[download] {n}/{len(pending)} {rec['requested']} {rec['download']} {rec['error']}", file=sys.stderr)
    os.replace(sfile + ".new", sfile)
    nfail = sum(1 for _ in csv.DictReader(open(sfile), delimiter="\t") if _["download"] != "ok")
    print(f"[download] shard {a.shard} finished, {nfail} failed (re-run to retry)", file=sys.stderr)
    sys.exit(0)   # failures are recorded, not fatal: re-submit to retry only those

# ---------------------------------------------------------------- merge
def cmd_merge(a):
    res = list(csv.DictReader(open(os.path.join(a.outdir, "resolution.tsv")), delimiter="\t"))
    dl = {}
    sd = os.path.join(a.outdir, "status")
    for f in sorted(os.listdir(sd)) if os.path.isdir(sd) else []:
        if f.endswith(".tsv"):
            for r in csv.DictReader(open(os.path.join(sd, f)), delimiter="\t"):
                dl[r["requested"]] = r
    cols = RES_COLS + ["download", "download_error", "fna_path", "genome_sha256"]
    from collections import Counter
    cnt = Counter()
    with open(os.path.join(a.outdir, "manifest.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t"); w.writeheader()
        for r in res:
            d = dl.get(r["requested"])
            if r["status"] in ("ok", "replaced"):
                state = d["download"] if d else "pending"
            else:
                state = "skipped"
            r.update(download=state, download_error=(d or {}).get("error", ""),
                     fna_path=os.path.abspath(os.path.join(a.outdir, "genomes", r["requested"] + ".fna.gz")) if state == "ok" else "",
                     genome_sha256=(d or {}).get("genome_sha256", ""))
            w.writerow(r); cnt[(r["status"], state)] += 1
    with open(os.path.join(a.outdir, "summary.txt"), "w") as fh:
        for k, v in sorted(cnt.items()):
            fh.write(f"{k[0]}\t{k[1]}\t{v}\n")
    print(open(os.path.join(a.outdir, "summary.txt")).read())
    unfinished = sum(v for (s, d), v in cnt.items() if d in ("failed", "pending"))
    sys.exit(1 if unfinished else 0)

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("resolve"); r.add_argument("accessions"); r.add_argument("outdir")
    r.add_argument("--replacement", choices=["any", "none"], default="any",
                   help="any: substitute retired accessions with their live successor/counterpart (flagged in manifest); none: just report them")
    d = sp.add_parser("download"); d.add_argument("outdir")
    d.add_argument("--include", default="genome", help="comma list of: " + ",".join(INCLUDE_MAP))
    d.add_argument("--shard", type=int, default=int(os.environ.get("SLURM_ARRAY_TASK_ID", 0)))
    d.add_argument("--nshards", type=int, default=1)
    d.add_argument("--threads", type=int, default=2)
    d.add_argument("--rps", type=float, default=3.0, help="max requests/s for THIS process (keep sum over shards <=5, or <=10 with NCBI_API_KEY)")
    m = sp.add_parser("merge"); m.add_argument("outdir")
    a = p.parse_args()
    {"resolve": cmd_resolve, "download": cmd_download, "merge": cmd_merge}[a.cmd](a)

if __name__ == "__main__":
    main()
