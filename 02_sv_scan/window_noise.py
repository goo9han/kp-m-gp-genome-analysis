#!/usr/bin/env python3
"""Empirical noise of the normalized depth of a 9-kb window (the size of the dip interior, 84-93 kb) genome-wide,
for a BAM: tile contigs >= 100 kb of the reference into 9-kb windows, samtools bedcov -Q 20 -> mean depth / genome median
(median from profile_<label>.json). Excludes windows within 20 kb of contig ends.
usage: window_noise.py <label> <bam> <ref.fa>"""
import sys, os, subprocess, json
import numpy as np
label, bam, ref = sys.argv[1:4]
SAM = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/samtools")
med = json.load(open(f"profile_{label}.json"))["median"]
bed = f"work_{label}_9kb.bed"
with open(bed, "w") as fh:
    for l in open(ref + ".fai"):
        c, L = l.split("\t")[0], int(l.split("\t")[1])
        if L < 100000: continue
        for s in range(20000, L - 20000 - 9000, 9000): fh.write(f"{c}\t{s}\t{s + 9000}\n")
out = subprocess.run([SAM, "bedcov", "-Q", "20", bed, bam], capture_output=True, text=True).stdout
v = np.array([int(l.split("\t")[3]) / 9000 / med for l in out.splitlines()])
q = np.percentile(v, [1, 2.5, 5, 25, 50, 75, 95, 97.5, 99])
res = dict(label=label, n=len(v), quantiles=dict(zip(["1", "2.5", "5", "25", "50", "75", "95", "97.5", "99"], q.tolist())),
           sd_central=float(np.std(v[(v > q[1]) & (v < q[7])])), mad_sd=float(1.4826 * np.median(np.abs(v - np.median(v)))))
json.dump(res, open(f"noise_{label}.json", "w")); print(res)
