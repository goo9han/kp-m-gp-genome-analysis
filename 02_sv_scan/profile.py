#!/usr/bin/env python3
"""100-bp bin depth along one contig, normalized by the genome median of the same BAM.
usage: profile.py <label> <bam> <ref.fa> <contig> [mapq]
genome median: samtools depth -a -Q <mapq> over contigs >= 20 kb of <ref.fa>, every 100th position (as win_het.py).
Output: profile_<label>.json {median, contig, length, bins:[start, mean_depth, norm]} (bins 1-based starts)."""
import sys, os, subprocess, json
import numpy as np
label, bam, ref, ctg = sys.argv[1:5]; Q = sys.argv[5] if len(sys.argv) > 5 else "20"
SAM = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/samtools")
fai = [(l.split("\t")[0], int(l.split("\t")[1])) for l in open(ref + ".fai")]
bed = f"work_{label}.bed"
open(bed, "w").write("".join(f"{c}\t0\t{L}\n" for c, L in fai if L >= 20000))
p = subprocess.Popen([SAM, "depth", "-a", "-Q", Q, "-b", bed, bam], stdout=subprocess.PIPE, text=True)
vals = [int(l.rsplit("\t", 1)[1]) for i, l in enumerate(p.stdout) if i % 100 == 0]; p.wait()
med = float(np.median(vals))
L = dict(fai)[ctg]
d = subprocess.run([SAM, "depth", "-a", "-Q", Q, "-r", ctg, bam], capture_output=True, text=True).stdout
dep = np.zeros(L + 1)
for l in d.splitlines():
    _, pos, v = l.split("\t"); dep[int(pos)] = int(v)
bins = []
for s in range(1, L + 1, 100):
    m = float(dep[s:min(L, s + 99) + 1].mean()); bins.append([s, m, m / med])
np.save(f"perbase_{label}.npy", dep)
json.dump(dict(label=label, bam=bam, ref=ref, contig=ctg, length=L, mapq=int(Q), median=med, n_median_positions=len(vals), bins=bins),
          open(f"profile_{label}.json", "w"))
print(label, "median", med, "contig", ctg, L)
