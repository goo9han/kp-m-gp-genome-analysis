#!/usr/bin/env python3
"""Per-site A/B haplotype bases on kp1 coordinates (mapping (b)).
usage: call_ab.py <label> <bam>
A = major allele, B = minor allele at heterozygous sites (same het rule as win_het.py) with MAF <= 0.45
(2:1 dosage: the single-copy allele is B). Homozygous sites (top allele fraction >= 0.9): A = B = top allele.
Everything else (uncallable, MAF > 0.45, 0.1-0.2 minor fraction, deletions '*') -> N in both.
This is dosage-based assignment, not read-backed phasing; check_linkage.py tests it with reads.
Output: ab_<label>.npz with per-scaffold uint8 arrays (0=N,1=A,2=C,3=G,4=T) A_<ctg>, B_<ctg>, and float16 maf_<ctg>."""
import os, subprocess, sys
from multiprocessing import Pool
import numpy as np
from scipy import stats

label, bam = sys.argv[1], sys.argv[2]
REF = "/Users/geongoohan/SynologyDrive/Brown/006_Fungal_WGS/kp_wgs_pangenome/qpcr_primer_design/validation/genomes/kp1_combined.fa"
BCF = os.path.expanduser("~/miniconda3/envs/bcftools_env/bin/bcftools")
SAM = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/samtools")
D_MIN, ERR, P = 14, 0.01, 1e-3
CODE = {"A": 1, "C": 2, "G": 3, "T": 4}
contigs = [(l.split("\t")[0], int(l.split("\t")[1])) for l in open(REF + ".fai")]
contigs = [(c, L) for c, L in contigs if L >= 20000]
_cache = {}


def is_het(b, n):
    k = (b, n)
    if k not in _cache:
        maf = b / n
        _cache[k] = (0.2 <= maf <= 0.8 and stats.binomtest(b, n, ERR, alternative="greater").pvalue < P
                     and stats.binomtest(n - b, n, 1 - ERR, alternative="less").pvalue < P)
    return _cache[k]


def run(ctg):
    c, L = ctg
    cmd = f"{BCF} mpileup -B -q 20 -Q 20 -I -a AD --max-depth 100000 -f {REF} -r {c} {bam} 2>/dev/null | " \
          f"{BCF} query -f '%POS\\t%REF\\t%ALT\\t[%AD]\\n'"
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout
    A = np.zeros(L + 1, np.uint8); B = np.zeros(L + 1, np.uint8); M = np.zeros(L + 1, np.float16)
    for line in out.splitlines():
        pos, rb, alt, ad = line.split("\t")
        ad = [int(x) for x in ad.split(",") if x != "."]
        n = sum(ad)
        if n < D_MIN or n > UMAX:
            continue
        pos = int(pos)
        alleles = [rb.upper()] + alt.split(",")
        order = sorted(range(len(ad)), key=lambda i: -ad[i])
        a = ad[order[0]]; top = alleles[order[0]]
        if top not in CODE:
            continue
        if a / n >= 0.9:
            A[pos] = B[pos] = CODE[top]; continue
        if len(ad) < 2:
            continue
        b = ad[order[1]]; sec = alleles[order[1]]
        if sec in CODE and b >= 2 and is_het(b, n) and b / n <= 0.45:
            A[pos] = CODE[top]; B[pos] = CODE[sec]; M[pos] = b / n
    return c, A, B, M


if __name__ == "__main__":
    import multiprocessing as mp
    mp.set_start_method("fork")
    MED = float(__import__("json").load(open(f"../win_{label}.json"))["median_depth"]); UMAX = 2.5 * MED
    globals().update(UMAX=UMAX)
    with Pool(6) as pool:
        res = pool.map(run, sorted(contigs, key=lambda x: -x[1]))
    d = {}
    for c, A, B, M in res:
        d[f"A_{c}"] = A; d[f"B_{c}"] = B; d[f"maf_{c}"] = M
    np.savez_compressed(f"ab_{label}.npz", **d)
    na = sum(int((x > 0).sum()) for k, x in d.items() if k.startswith("A_"))
    nh = sum(int(((d[k] > 0) & (d[k] != d["B_" + k[2:]])).sum()) for k in d if k.startswith("A_"))
    print(label, "called", na, "A!=B", nh, flush=True)
