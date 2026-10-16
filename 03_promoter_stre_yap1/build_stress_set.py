#!/usr/bin/env python3
"""STRESS promoter set = organic-acid true-ATG promoters (60 loci) + H2O2 true-ATG promoters (14 loci),
de-duplicated by kp1 funannotate gene (organic-acid copy kept). Loci with a FAILED strain are dropped.
Writes work/stress.fa (gid__strain), work/stress_meta.json (labels, contig ranges, flags) and
work/stress_loci.json (all kp1 gene IDs touched, used to exclude them from the BACKGROUND draw).
"""
import csv, json, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
STRAINS = ["kp1", "kpmtx", "C89-2", "ID1LA", "IDNYC", "QC1", "W26G", "kpgp3", "kpgp6"]


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


def parse_range(r):
    m = re.match(r"(.+):(\d+)-(\d+)\(([+-])\)$", r)
    return m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)


meta, fa = {}, []
# ---- organic acid
X = json.load(open(f"{ROOT}/organic_acid_analysis/true_atg_rescreen/work/extraction.json"))
oa = read_fasta(f"{ROOT}/organic_acid_analysis/true_atg_rescreen/OA_promoters_trueATG_9strains.fa")
rows = {}
for r in X["rows"]:
    rows.setdefault(r["kp1_gene"], {})[r["strain"]] = r
for gid, rs in rows.items():
    if any(rs[s].get("extraction") == "FAILED" for s in STRAINS):
        meta[gid] = dict(label=X["locus_label"][gid], source="organic_acid", used=False,
                         reason="FAILED strain(s): " + ",".join(s for s in STRAINS if rs[s].get("extraction") == "FAILED"))
        continue
    flags = [f"{s}:{rs[s]['extraction']}" for s in STRAINS if rs[s]["extraction"] != "FULL"]
    meta[gid] = dict(label=X["locus_label"][gid], source="organic_acid", used=True, flags=flags,
                     ranges={s: parse_range(rs[s]["promoter_range"]) for s in STRAINS})
    for s in STRAINS:
        fa.append(f">{gid}__{s}\n{oa[f'{gid}__{s}']}")
# ---- H2O2
src = open(f"{ROOT}/h2o2_analysis/true_atg_all/h2o2_true_atg_extract.py").read()
lab2gid = dict(re.findall(r'\("(\w+)", "\w+", \[[^\]]*\], "(FUN_\d+)"', src))
h2 = read_fasta(f"{ROOT}/h2o2_analysis/true_atg_all/H2O2_promoters_trueATG_9strains.fa")
hlog = {}
for r in csv.DictReader(open(f"{ROOT}/h2o2_analysis/true_atg_all/true_atg_extraction_log.tsv"), delimiter="\t"):
    hlog.setdefault(r["locus"], {})[r["strain"]] = r
for lab, gid in lab2gid.items():
    if gid in meta:
        meta[gid].setdefault("also_in", []).append(f"H2O2:{lab}")
        continue
    rs = hlog[lab]
    flags = [f"{s}:{rs[s]['note']}" for s in STRAINS if "padded" in rs[s]["note"]]
    meta[gid] = dict(label=lab, source="h2o2", used=True, flags=flags,
                     ranges={s: parse_range(rs[s]["new_promoter_contig_range"]) for s in STRAINS})
    for s in STRAINS:
        fa.append(f">{gid}__{s}\n{h2[f'{lab}__{s}']}")
os.makedirs(f"{HERE}/work", exist_ok=True)
open(f"{HERE}/work/stress.fa", "w").write("\n".join(fa) + "\n")
json.dump(meta, open(f"{HERE}/work/stress_meta.json", "w"), indent=1)
json.dump(sorted(meta), open(f"{HERE}/work/stress_loci.json", "w"))
used = [g for g, m in meta.items() if m["used"]]
print(f"STRESS loci: {len(meta)} total, {len(used)} used; dropped: " +
      "; ".join(f"{m['label']} ({m['reason']})" for m in meta.values() if not m["used"]))
print("shared OA/H2O2:", {g: m["also_in"] for g, m in meta.items() if "also_in" in m})
