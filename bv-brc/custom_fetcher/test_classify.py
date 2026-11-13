import sys; sys.path.insert(0, ".")
import fetch_genomes as f
live=lambda a,s="current":{"accession":a,"assembly_info":{"assembly_status":s}}
rep={"GCF_000000001.1":live("GCF_000000001.1","suppressed"),"GCA_000000001.1":live("GCA_000000001.1")}
rl=lambda a:[{"genbank_accession":"GCA_000000001.1","refseq_accession":"GCF_000000001.1","identical":True}]
r=f.classify("GCF_000000001.1",rep,rl); assert (r["status"],r["resolved"],r["resolution"])==("replaced","GCA_000000001.1","other_database_counterpart"); print(r["status"],r["resolved"],r["resolution"],"|",r["note"])
r=f.classify("GCF_000000001.1",{"GCF_000000001.1":live("GCF_000000001.1","suppressed")},lambda a:[]); assert r["status"]=="retired"; print(r["status"],"|",r["note"])
r=f.classify("GCF_000000001.1",{},rl,allow_replace=False); assert r["status"]=="retired"; print(r["status"],"|",r["note"])
print("unit tests passed")
# retry layer survives proxy-style protocol errors (BadStatusLine) and then succeeds
import http.client, urllib.request, io
calls = {"n": 0}
class R(io.BytesIO):
    status = 200
    def __enter__(self): return self
    def __exit__(self, *a): pass
def fake(req, timeout=None):
    calls["n"] += 1
    if calls["n"] < 3: raise http.client.BadStatusLine("CONNECT x HTTP/1.0")
    return R(b'{"ok": 1}')
f.time.sleep = lambda s: None; f.LIMITER.wait = lambda: None
urllib.request.urlopen = fake
assert f.http_get("https://x/y") == (200, {"ok": 1}) and calls["n"] == 3
print("retry test passed")
