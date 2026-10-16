#!/usr/bin/env python3
"""Pull the Kp-confident nuclear reads out of a mapped BAM and write them as single-end FASTQ.

Reference is the deposited Kp-GP genome, nuclear and mitochondrial sequence together. The nuclear
assembly is its own file, so "nuclear" is a file boundary here rather than a name test, and the only
condition applied to the contig is the length one.

  primary alignment, mapping quality at least 20, aligned length at least 100 bp, identity at least
  97 per cent, on a contig of at least 20 kb of the deposited nuclear assembly.

Sequences are written as the BAM stores them, which is reference-oriented. The scoring path realigns
single-end, so orientation does not affect the score, and writing them this way is what reproduced the
earlier read sets byte for byte.

usage: extract_kp_confident.py <run> <bam> <out.fq>
"""
import os, re, sys, subprocess

KP = "/Users/geongoohan/SynologyDrive/Brown/006_Fungal_WGS/kp_wgs_pangenome"
SAM = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/samtools")
NUCFAI = f"{KP}/public_survey/rigp_ref/nuclear_only.fai"
NUC = {c for c, L in (l.split("\t")[:2] for l in open(NUCFAI)) if int(L) >= 20000}


def run(label, bam, out):
    p = subprocess.Popen([SAM, "view", "-F", "0x904", "-q", "20", bam],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    n = 0
    with open(out, "w") as fh:
        for l in p.stdout:
            f = l.split("\t", 12)
            if f[2] not in NUC:
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
            if seq == "*" or len(seq) != len(qual):
                continue
            fh.write(f"@{label}_{n}\n{seq}\n+\n{qual}\n")
            n += 1
    p.wait()
    print(f"  {label:<12} Kp-confident nuclear reads {n:>9,}  ->  {out}", flush=True)
    return n


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], sys.argv[3])
