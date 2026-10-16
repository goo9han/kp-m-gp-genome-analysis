#!/usr/bin/env python3
"""1-kb bin depth over the whole assembly for one existing BAM (no remapping).
usage: depth_bins.py <label> <bam> <ref.fa> <median|auto>
- Bins: contigs >= 20 kb of <ref.fa> (mitochondrial contig excluded), tiled from position 0 in 1-kb steps; the last
  partial bin is kept with its true width.
- Depth: samtools bedcov -Q 20 -j (MAPQ >= 20, deletions/ref-skips not counted = same counting as samtools depth used by
  ahp1_cnv/profile.py) -> mean depth per bin.
- Normalization: the strain's genome median exactly as ahp1_cnv/profile.py (samtools depth -a -Q 20 over contigs >= 20 kb,
  every 100th position). <median> is copied from the existing ahp1_cnv/profile_<label>.json; 'auto' recomputes it that way
  (only for BAMs without a stored profile, i.e. the GP-on-kp1 BAMs).
Output: work/bins_<label>.tsv (contig, start0, end, mean_depth, norm) and work/median_<label>.txt"""
import sys, os, subprocess
import numpy as np
label, bam, ref, med = sys.argv[1:5]
SAM = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/samtools")
HERE = os.path.dirname(os.path.abspath(__file__)); W = f"{HERE}/work"; os.makedirs(W, exist_ok=True)
fai = [(l.split("\t")[0], int(l.split("\t")[1])) for l in open(ref + ".fai")]
ctgs = [(c, L) for c, L in fai if L >= 20000 and not c.endswith("_MT")]   # Kp_GP_MT / Kp_M_MT
if med == "auto":
    bed = f"{W}/ge20k_{label}.bed"
    open(bed, "w").write("".join(f"{c}\t0\t{L}\n" for c, L in fai if L >= 20000))   # same set as profile.py (mt included there)
    p = subprocess.Popen([SAM, "depth", "-a", "-Q", "20", "-b", bed, bam], stdout=subprocess.PIPE, text=True)
    vals = [int(l.rsplit("\t", 1)[1]) for i, l in enumerate(p.stdout) if i % 100 == 0]; p.wait()
    med = float(np.median(vals))
else:
    med = float(med)
open(f"{W}/median_{label}.txt", "w").write(f"{med}\n")
bed = f"{W}/bins1k_{os.path.basename(ref)}.bed"
if not os.path.exists(bed):
    with open(bed, "w") as fh:
        for c, L in ctgs:
            for s in range(0, L, 1000): fh.write(f"{c}\t{s}\t{min(L, s + 1000)}\n")
out = subprocess.run([SAM, "bedcov", "-Q", "20", "-j", bed, bam], capture_output=True, text=True, check=True).stdout
with open(f"{W}/bins_{label}.tsv", "w") as fh:
    for l in out.splitlines():
        c, s, e, tot = l.split("\t"); s, e = int(s), int(e); m = int(tot) / (e - s)
        fh.write(f"{c}\t{s}\t{e}\t{m:.3f}\t{m / med:.5f}\n")
print(label, "median", med, "bins", len(out.splitlines()))
