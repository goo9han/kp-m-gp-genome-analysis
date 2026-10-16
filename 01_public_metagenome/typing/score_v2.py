#!/usr/bin/env python3
"""Score a read set at the informative positions of build_sites_twoisolate.py and keep enough detail
for the mixture controls and for the genome-wide uniformity check.

Per read it records how many observations carry the Kp-M base and how many carry the Kp-GP allele, so a
mixture draw is an index-and-sum rather than a rescan. Per contig it records the same totals.

Filters are the published ones: bwa mem to kp1_combined.fa, primary alignments only, mapping quality at
least 20, aligned length at least 100 bp, identity at least 97 per cent, and only bases of quality 20 or
higher.

usage: score_detail.py <label...>      reads <label>.fq from the known work directories
Writes detail_<label>.npz with per-read arrays m, g and per-contig arrays contig, cm, cg.
"""
import os, re, sys, subprocess, collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
KP = os.path.dirname(os.path.dirname(HERE))
SAM = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/samtools")
BWA = os.path.expanduser("~/miniconda3/envs/pangenome/bin/bwa")
K1 = f"{KP}/qpcr_primer_design/validation/genomes/kp1_combined.fa"
FQDIRS = [HERE + "/work"]
S = {k: v for k, v in np.load(f"{KP}/public_survey/nine_panel_test/sites_twoisolate.npz").items()}


def fq_for(label):
    for d in FQDIRS:
        p = f"{d}/{label}.fq"
        if os.path.exists(p):
            return p
    raise SystemExit(f"no FASTQ for {label}")


def score(label):
    fq = fq_for(label)
    p = subprocess.Popen(f"{BWA} mem -t 8 {K1} {fq} 2>/dev/null | {SAM} view -F 0x904 -q 20 -",
                         shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    per = collections.defaultdict(lambda: [0, 0])          # read name -> [M obs, GP obs]
    byc = collections.defaultdict(lambda: [0, 0])          # contig -> [M obs, GP obs]
    for l in p.stdout:
        f = l.split("\t", 12)
        c = f[2]
        if f"inf_{c}" not in S:
            continue
        ops = re.findall(r"(\d+)([MIDNSHP=X])", f[5])
        aln = sum(int(a) for a, o in ops if o in "M=X")
        if aln < 100:
            continue
        nm = None
        for t in f[11].split("\t"):
            if t.startswith("NM:i:"):
                nm = int(t[5:]); break
        if nm is None or 1 - nm / aln < 0.97:
            continue
        seq, qual = f[9], f[10]
        rp, qp = int(f[3]) - 1, 0
        inf = S[f"inf_{c}"]; gp = S[f"gp_{c}"]; rf = S[f"ref_{c}"]
        for a, o in ops:
            a = int(a)
            if o in "M=X":
                rr = np.arange(rp, rp + a)
                msk = inf[rr]
                if msk.any():
                    qq = np.arange(qp, qp + a)
                    for r0, q0 in zip(rr[msk], qq[msk]):
                        if ord(qual[q0]) - 33 < 20:
                            continue
                        b = "NACGT".index(seq[q0]) if seq[q0] in "ACGT" else 0
                        if b == rf[r0]:
                            per[f[0]][0] += 1; byc[c][0] += 1
                        elif b == gp[r0]:
                            per[f[0]][1] += 1; byc[c][1] += 1
                rp += a; qp += a
            elif o in "DN":
                rp += a
            elif o in "IS":
                qp += a
    p.wait()
    names = sorted(per)
    m = np.array([per[k][0] for k in names], np.int32)
    g = np.array([per[k][1] for k in names], np.int32)
    cs = sorted(byc)
    np.savez_compressed(f"{HERE}/detail_{label}.npz", m=m, g=g,
                        contig=np.array(cs), cm=np.array([byc[c][0] for c in cs], np.int64),
                        cg=np.array([byc[c][1] for c in cs], np.int64))
    tm, tg = int(m.sum()), int(g.sum())
    print(f"  {label:<12} reads with obs {len(names):>8,}   M {tm:>9,}  GP {tg:>9,}   "
          f"GP {100 * tg / (tm + tg):.2f}%", flush=True)


if __name__ == "__main__":
    for lab in sys.argv[1:]:
        score(lab)
