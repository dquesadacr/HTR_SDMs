#!/usr/bin/env python3
"""01_pull_gbif.py — broad-scale occurrences for the truncation PoC. Run on PIK login node (internet)."""
import json, time, csv, sys, urllib.request

SP = {"Arnica montana": 5405976, "Cicerbita alpina": 3087634, "Meum athamanticum": 6026490}
MAXREC = 9000   # GBIF REST offset+limit cap is 10000
BASE = "https://api.gbif.org/v1/occurrence/search"
Q = "&" + urllib.parse.quote(
    "hasCoordinate=true&occurrenceStatus=PRESENT&decimalLatitude=34,62"
    "&decimalLongitude=-12,25&coordinateUncertaintyInMeters=0,20000"
    "&year=1970,2000", safe="=&:,-")

def fetch(url):
    for a in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as r: return json.load(r)
        except Exception as e:
            print("retry", a, e, file=sys.stderr); time.sleep(2 * (a + 1))
    raise RuntimeError("gbif fetch failed: " + url)

for name, key in SP.items():
    rows, off = {}, 0
    total = None
    while total is None or off < min(total, MAXREC):
        d = fetch(f"{BASE}?taxonKey={key}{Q}&limit=1000&offset={off}")
        if total is None:
            total = d.get("count") or d.get("matchedRecords", 0)
            print(f"{name}: matchedRecords={total}")
        recs = d.get("results", [])
        if not recs: break
        for r in recs:
            try: k = (round(r["decimalLatitude"], 2), round(r["decimalLongitude"], 2))
            except Exception: continue
            rows[k] = (r["decimalLongitude"], r["decimalLatitude"], int(r.get("year", 0)))
        off += 1000
        time.sleep(0.2)
    out = f"gbif_{name.replace(' ', '_')}.csv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["lon", "lat", "year"])
        for v in sorted(rows.values()): w.writerow(v)
    print(name, "->", len(rows), "deduped rows ->", out)
print("done")
