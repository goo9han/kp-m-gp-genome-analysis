#!/usr/bin/env python3
"""All copies of the 316-bp HOG1-AHP1 boundary repeat in the kpgp3 (GP) and kp1 (mouse) assemblies, and every pair of copies on
the same contig whose inner gap is 5-30 kb. For each pair: orientation (direct = same strand, the configuration that can
delete the intervening sequence by recombination; inverted = opposite strands), identity between the two copies, genes in the
gap, and GP / mouse depth over the gap (1-kb bins from depth_bins.py, normalized to own genome median).
Query: kpgp3 scaffold_18:79,562-79,877 (left copy). blastn -dust no, e <= 1e-10, alignment >= 200 bp counts as a copy.
Output: repeat_pairs.tsv, repeats/copies_<asm>.tsv"""
import os, subprocess, re, collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.dirname(HERE); W = f"{HERE}/work"
B = os.path.expanduser("~/miniconda3/envs/buscoenv/bin")
G = f"{R}/qpcr_primer_design/validation/genomes"
ASM = {"kpgp3": (f"{G}/kpgp3_combined.fa", f"{G}/kpgp3_db",
                 "/Users/geongoohan/SynologyDrive/Brown/006_Fungal_WGS/kpgp3/results/funannotate/Kazachstania_pintolopesii.gff3",
                 "/Users/geongoohan/SynologyDrive/Brown/006_Fungal_WGS/kpgp3/results/funannotate/Kazachstania_pintolopesii.emapper.annotations"),
       "kp1": (f"{G}/kp1_combined.fa", f"{G}/kp1_db",
               "/Users/geongoohan/SynologyDrive/Brown/006_Fungal_WGS/kp1/scaffolding_medaka_racon_pilon/funannotate/Kazachstania_pintolopesii.gff3",
               "/Users/geongoohan/SynologyDrive/Brown/006_Fungal_WGS/kp1/scaffolding_medaka_racon_pilon/funannotate/Kazachstania_pintolopesii.emapper.annotations")}
DEPTH = {"kpgp3": (["kpgp3_self", "kpgp6_on3"], ["kp1_on3", "kpmtx_on3", "C89-2_on3", "ID1LA_on3", "IDNYC_on3", "QC1_on3", "W26G_on3"]),
         "kp1": (["kpgp3_on_k1", "kpgp6_on_k1"], ["kp1_self_k1"])}
Q = f"{HERE}/repeats/rep316_copyL.fa"
HOG1_AHP1 = ("scaffold_18", 79562, 93667)


def genes_in(gff, emp, c, s, e):
    if not gff: return "annotation not available"
    ann = {}
    if emp and os.path.exists(emp):
        for l in open(emp):
            if l.startswith("#"): continue
            f = l.rstrip("\n").split("\t"); ann[f[0].split("-T")[0]] = f[8] or f[7][:40]
    out = []
    for l in open(gff):
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] != c or f[2] != "gene": continue
        gs, ge = int(f[3]), int(f[4])
        if ge < s or gs > e: continue
        a = dict(x.split("=", 1) for x in f[8].strip(";").split(";") if "=" in x)
        nm = a.get("Name", "") or ann.get(a["ID"], "")
        out.append(f"{a['ID']}{'(' + nm + ')' if nm else ''}")
    return ",".join(out) if out else "none"


def depth_between(asm, c, s, e):
    gp, mo = DEPTH[asm]; res = {}
    for lab in gp + mo:
        p = f"{W}/bins_{lab}.tsv"
        v = [float(x.split("\t")[4]) for x in open(p) if x.split("\t")[0] == c and int(x.split("\t")[1]) >= s and int(x.split("\t")[2]) <= e]
        res[lab] = float(np.mean(v)) if v else float("nan")
    return res, gp, mo


def identity(fa, c1, c2):
    f1, f2 = f"{W}/_a.fa", f"{W}/_b.fa"
    for (c, s, e), f in ((c1, f1), (c2, f2)):
        open(f, "w").write(subprocess.run([f"{B}/samtools", "faidx", fa, f"{c}:{s}-{e}"], capture_output=True, text=True).stdout)
    o = subprocess.run([f"{B}/blastn", "-query", f1, "-subject", f2, "-dust", "no", "-outfmt", "6 pident length"],
                       capture_output=True, text=True).stdout.splitlines()
    return max((float(x.split()[0]), int(x.split()[1])) for x in o) if o else (0.0, 0)


rows = []
for asm, (fa, db, gff, emp) in ASM.items():
    raw = subprocess.run([f"{B}/blastn", "-query", Q, "-db", db, "-dust", "no", "-evalue", "1e-10", "-max_target_seqs", "5000",
                          "-max_hsps", "5000", "-outfmt", "6 sseqid pident length sstart send evalue"], capture_output=True, text=True).stdout
    cp = []
    for l in raw.splitlines():
        c, pid, ln, a, b, ev = l.split("\t"); a, b = int(a), int(b)
        if int(ln) >= 200: cp.append((c, min(a, b), max(a, b), "+" if b > a else "-", float(pid)))
    cp.sort()
    with open(f"{HERE}/repeats/copies_{asm}.tsv", "w") as fh:
        fh.write("contig\tstart\tend\tstrand_vs_query\tpident_vs_query\n")
        for x in cp: fh.write("\t".join(map(str, x)) + "\n")
    by = collections.defaultdict(list)
    for x in cp: by[x[0]].append(x)
    n_pairs = 0
    for c, L in by.items():
        for i in range(len(L)):
            for j in range(i + 1, len(L)):
                a, b = L[i], L[j]; gap = b[1] - a[2] - 1
                if not 5000 <= gap <= 30000: continue
                n_pairs += 1
                pid, aln = identity(fa, (c, a[1], a[2]), (c, b[1], b[2]))
                d, gp, mo = depth_between(asm, c, a[2], b[1])
                is_ctrl = asm == "kpgp3" and c == HOG1_AHP1[0] and a[1] <= HOG1_AHP1[1] + 5 and b[2] >= HOG1_AHP1[2] - 5
                rows.append(dict(assembly=asm, contig=c, copy1=f"{a[1]}-{a[2]}", copy1_strand=a[3], copy2=f"{b[1]}-{b[2]}",
                                 copy2_strand=b[3], orientation="direct" if a[3] == b[3] else "inverted", gap_bp=gap,
                                 copy_pair_identity=f"{pid:.1f} ({aln} bp)",
                                 gp_depth_gap=" / ".join(f"{d[g]:.2f}" for g in gp), mouse_depth_gap=" / ".join(f"{d[m]:.2f}" for m in mo),
                                 gp_min=min(d[g] for g in gp), genes_in_gap=genes_in(gff, emp, c, a[2], b[1]),
                                 note="HOG1-AHP1 pair (positive control)" if is_ctrl else ""))
    print(asm, "copies", len(cp), "pairs 5-30 kb", n_pairs)
cols = ["assembly", "contig", "copy1", "copy1_strand", "copy2", "copy2_strand", "orientation", "gap_bp", "copy_pair_identity",
        "gp_depth_gap", "mouse_depth_gap", "genes_in_gap", "note"]
with open(f"{HERE}/repeat_pairs.tsv", "w") as fh:
    fh.write("# Pairs of copies of the 316-bp HOG1-AHP1 boundary repeat (query kpgp3 scaffold_18:79,562-79,877; blastn -dust no, e<=1e-10, "
             "aln>=200 bp) on the same contig with a 5-30 kb inner gap. Depth = mean normalized 1-kb bin depth over the gap "
             "(kpgp3 assembly: GP = Kp-GP / Kp-GP2, mouse = Kp-M / Kp-MTX / C89-2 / ID1LA / IDNYC / QC1 / W26G; "
             "kp1 assembly: GP reads on kp1, mouse = Kp-M self).\n")
    fh.write("\t".join(cols) + "\n")
    for r in sorted(rows, key=lambda r: (r["assembly"] != "kpgp3", r["gp_min"])):
        fh.write("\t".join(str(r[k]) for k in cols) + "\n")
print(open(f"{HERE}/repeat_pairs.tsv").read())
