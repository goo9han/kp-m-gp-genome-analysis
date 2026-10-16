#!/usr/bin/env python3
"""BACKGROUND promoter set: 300 random kp1 funannotate genes (seed 20260925) not in the STRESS set,
extracted with the true-ATG rule of organic_acid_analysis/true_atg_rescreen/oa_true_atg_extract.py:
kp1 = funannotate CDS start; other strains = miniprot rank-1 hit starting at kp1 residue 1 with
identity >= 0.90; promoter = genome -700..-1 on the coding strand; spliced ORF without a stop in the
first 100 codons (or a terminal stop for shorter proteins). Loci with any failed strain are dropped.
"""
import json, os, random, re, subprocess
from collections import defaultdict
from Bio.Seq import Seq

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GENOME_DIR = f"{ROOT}/rnaseq_prep/ani/nuclear_clean"
GTF = f"{ROOT}/rnaseq_prep/reference/kp1_funannotate.gtf"
MINIPROT = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/miniprot")
STRAINS = ["kp1", "kpmtx", "C89-2", "ID1LA", "IDNYC", "QC1", "W26G", "kpgp3", "kpgp6"]
PROM_LEN, MIN_CODONS, N_BG, SEED = 700, 100, 300, 20260925
rc = lambda s: s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


translate = lambda nt: str(Seq(nt[:len(nt) // 3 * 3]).translate())
genomes = {s: read_fasta(f"{GENOME_DIR}/{s}_nuclear_clean.fa") for s in STRAINS}
cds = defaultdict(list)
for line in open(GTF):
    f = line.rstrip("\n").split("\t")
    if len(f) > 8 and f[2] == "CDS":
        cds[re.search(r'gene_id "([^"]+)"', f[8]).group(1)].append((f[0], int(f[3]), int(f[4]), f[6]))

stress = set(json.load(open(f"{HERE}/work/stress_loci.json")))
pool = sorted(g for g, segs in cds.items() if g not in stress and segs[0][0] in genomes["kp1"])
pick = random.Random(SEED).sample(pool, N_BG)


def spliced(g, segs, strand):
    segs = sorted(segs, key=lambda x: x[1])
    nt = "".join(g[segs[0][0]][a - 1:b] for _, a, b, _ in segs)
    return nt if strand == "+" else rc(nt)


def orf_ok(nt):
    p = translate(nt)
    i = p.find("*")
    return i == -1 or i >= MIN_CODONS or (i == len(p) - 1)


os.makedirs(f"{HERE}/work/bg", exist_ok=True)
prot, kp1 = {}, {}
for gid in pick:
    segs = cds[gid]
    contig, strand = segs[0][0], segs[0][3]
    nt = spliced(genomes["kp1"], segs, strand)
    p = translate(nt).rstrip("*")
    if "*" in p or len(p) < 30:
        continue
    seq = genomes["kp1"][contig] if strand == "+" else rc(genomes["kp1"][contig])
    k = (min(x[1] for x in segs) - 1) if strand == "+" else len(seq) - max(x[2] for x in segs)
    if seq[k:k + 3] != "ATG":
        continue
    prot[gid] = p
    kp1[gid] = dict(contig=contig, strand=strand, k=k, ok=orf_ok(nt))
with open(f"{HERE}/work/bg/kp1_bg_proteins.faa", "w") as fh:
    for g, p in prot.items():
        fh.write(f">{g}\n{p}\n")

mp = {}
for s in STRAINS[1:]:
    gff = f"{HERE}/work/bg/miniprot_{s}.gff"
    if not os.path.exists(gff):
        open(gff, "w").write(subprocess.run([MINIPROT, "-t", "8", "--gff", f"{GENOME_DIR}/{s}_nuclear_clean.fa",
                                             f"{HERE}/work/bg/kp1_bg_proteins.faa"], capture_output=True, text=True, check=True).stdout)
    mrna, cdsf = {}, defaultdict(list)
    for line in open(gff):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9:
            continue
        if f[2] == "mRNA" and "Rank=1;" in f[8]:
            t = re.search(r"Target=(\S+) (\d+) (\d+)", f[8])
            mrna[t.group(1)] = dict(id=re.search(r"ID=([^;]+)", f[8]).group(1), contig=f[0], start=int(f[3]), end=int(f[4]),
                                    strand=f[6], t0=int(t.group(2)), ident=float(re.search(r"Identity=([\d.]+)", f[8]).group(1)))
        elif f[2] == "CDS":
            cdsf[re.search(r"Parent=([^;]+)", f[8]).group(1)].append((f[0], int(f[3]), int(f[4]), f[6]))
    for m in mrna.values():
        m["cds"] = cdsf[m["id"]]
    mp[s] = mrna

kept, log, fa = [], [], []
for gid in prot:
    recs, fail = {}, None
    for s in STRAINS:
        if s == "kp1":
            contig, strand, k, ok = kp1[gid]["contig"], kp1[gid]["strand"], kp1[gid]["k"], kp1[gid]["ok"]
        else:
            m = mp[s].get(gid)
            if m is None or m["t0"] != 1 or m["ident"] < 0.90:
                fail = f"{s}: " + ("no hit" if m is None else f"Target start {m['t0']}, identity {m['ident']:.3f}"); break
            contig, strand = m["contig"], m["strand"]
            L = len(genomes[s][contig])
            k = m["start"] - 1 if strand == "+" else L - m["end"]
            ok = orf_ok(spliced(genomes[s], m["cds"], strand)) if m["cds"] else False
        seq = genomes[s][contig] if strand == "+" else rc(genomes[s][contig])
        if seq[k:k + 3] != "ATG":
            fail = f"{s}: start codon {seq[k:k + 3]}"; break
        if k < PROM_LEN:
            fail = f"{s}: contig ends {k} bp upstream"; break
        L = len(seq)
        rng = (contig, k - PROM_LEN + 1, k, "+") if strand == "+" else (contig, L - k + 1, L - k + PROM_LEN, "-")
        recs[s] = dict(prom=seq[k - PROM_LEN:k], range=rng, orf_ok=ok)
    if fail:
        log.append(f"{gid}\tDROPPED\t{fail}")
        continue
    kept.append(gid)
    log.append(f"{gid}\tKEPT\t" + ("" if all(r["orf_ok"] for r in recs.values()) else "ORF flag: " +
                                     ",".join(s for s in STRAINS if not recs[s]["orf_ok"])))
    for s in STRAINS:
        fa.append(f">{gid}__{s}\n{recs[s]['prom']}")
    json.dump({s: recs[s]["range"] for s in STRAINS}, open(f"{HERE}/work/bg/{gid}.range.json", "w"))
open(f"{HERE}/work/bg/BG_promoters_trueATG_9strains.fa", "w").write("\n".join(fa) + "\n")
open(f"{HERE}/work/bg/extraction_log.tsv", "w").write("gene\tstatus\tnote\n" + "\n".join(log) + "\n")
print(f"picked {N_BG}, protein ok {len(prot)}, kept {len(kept)}")
