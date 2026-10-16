#!/usr/bin/env python3
"""Check that the own PWM scorer reproduces FIMO's MSN2/MSN4 hit set (p<1e-3) exactly on the STRESS
promoters of all 9 strains (plan requirement), and on the BACKGROUND promoters."""
import os, subprocess, sys
from collections import defaultdict
from stre_core import *

FIMO = os.path.expanduser("~/miniconda3/envs/memeenv/bin/fimo")
MEME = f"{ROOT}/organic_acid_analysis/motifs_organic_acid.meme"
bg = load_bg()
M = read_meme(MEME)
pwms = [PWM(n, *M[n], bg) for n in ("MSN2_STRE", "MSN4_STRE")]
for p in pwms:
    print(f"{p.name}: score cutoff {p.cut:.4f} (p<1e-3)")
report = []
for tag, fa in (("STRESS", f"{HERE}/work/stress.fa"), ("BACKGROUND", f"{HERE}/work/bg/BG_promoters_trueATG_9strains.fa")):
    out = f"{HERE}/work/fimo_{tag}"
    if not os.path.exists(f"{out}/fimo.tsv"):
        subprocess.run([FIMO, "--thresh", "1e-3", "--max-stored-scores", "10000000", "--bfile",
                        f"{ROOT}/h2o2_analysis/meme_fimo/bg.txt", "--motif", "MSN2_STRE", "--motif", "MSN4_STRE",
                        "--oc", out, MEME, fa], check=True, capture_output=True)
    fimo = set()
    for line in open(f"{out}/fimo.tsv"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 10 or f[0] not in ("MSN2_STRE", "MSN4_STRE"):
            continue
        fimo.add((f[2], f[0], int(f[3]) - 1, int(f[4]) - 1, f[5]))
    own = set()
    for name, seq in read_fasta(fa).items():
        for p in pwms:
            for a, b, st, sc in p.hits(seq):
                own.add((name, p.name, a, b, st))
    only_f, only_o = fimo - own, own - fimo
    line = f"{tag}: FIMO hits {len(fimo)}, own hits {len(own)}, only FIMO {len(only_f)}, only own {len(only_o)}"
    print(line); report.append(line)
    for x in sorted(only_f)[:5] + sorted(only_o)[:5]:
        print("   ", x)
open(f"{HERE}/scorer_validation.txt", "w").write("\n".join(report) + "\n")
