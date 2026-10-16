#!/usr/bin/env python3
"""Step 4: every annotated feature (kpgp3 funannotate) overlapping the deleted interval and +-5 kb, with eggNOG
Sc ortholog text and GO-based category flags (direct + propagated terms as written in the GFF Ontology_term field).
Deleted interval = the sequence between the two copies of the flanking repeat that ONT deletion reads skip
(scaffold_18:79,562-93,351 or, equivalently, 79,878-93,667; 13,790 bp) - see ahp1_cnv_verdict.md."""
import re, csv
GFF = "/Users/geongoohan/SynologyDrive/Brown/006_Fungal_WGS/kpgp3/results/funannotate/Kazachstania_pintolopesii.gff3"
EMP = "/Users/geongoohan/SynologyDrive/Brown/006_Fungal_WGS/kpgp3/results/funannotate/Kazachstania_pintolopesii.emapper.annotations"
DEL = (79562, 93667)           # union of the two equivalent placements (one repeat copy is retained)
WIN = (DEL[0] - 5000, DEL[1] + 5000)
GO = {"stress (GO:0006950 response to stress)": "GO:0006950", "oxidative (GO:0006979 response to oxidative stress)": "GO:0006979",
      "pH (GO:0030003 / GO:0006885 / GO:0051453)": ("GO:0030003", "GO:0006885", "GO:0051453"),
      "transmembrane transport (GO:0055085)": "GO:0055085"}
emp = {}
for l in open(EMP):
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t"); emp[f[0].split("-T")[0]] = (f[7], f[8])
genes = {}
for l in open(GFF):
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t")
    if len(f) < 9 or f[0] != "scaffold_18": continue
    s, e = int(f[3]), int(f[4])
    if e < WIN[0] or s > WIN[1]: continue
    a = dict(x.split("=", 1) for x in f[8].strip(";").split(";") if "=" in x)
    if f[2] == "gene": genes[a["ID"]] = dict(gene_id=a["ID"], start=s, end=e, strand=f[6], name=a.get("Name", ""))
    elif f[2] in ("mRNA", "tRNA"):
        g = genes[a["Parent"]]; g["type"] = f[2]; g["product"] = a.get("product", ""); g["go"] = a.get("Ontology_term", "")
with open("ahp1_region_genes.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["gene_id", "type", "contig", "start", "end", "strand", "position_vs_deletion", "name", "product",
                "Sc_ortholog_eggNOG", "eggNOG_preferred_name"] + list(GO) + ["note"])
    for g in sorted(genes.values(), key=lambda g: g["start"]):
        pos = "INSIDE" if g["start"] >= DEL[0] and g["end"] <= DEL[1] else ("PARTIAL" if g["end"] >= DEL[0] and g["start"] <= DEL[1] else "outside (flank)")
        desc, pref = emp.get(g["gene_id"], ("", ""))
        flags = []
        for k, v in GO.items():
            v = v if isinstance(v, tuple) else (v,)
            flags.append("Y" if any(t in g.get("go", "") for t in v) else ("-" if g.get("type") == "mRNA" else "n/a"))
        note = ""
        if g["name"] == "HOG1": note = "eggNOG seed is not an Sc protein; Sc HOG1 = YLR113W by name"
        if g.get("type") == "tRNA": note = "tRNA gene; no eggNOG entry"
        w.writerow([g["gene_id"], g.get("type", ""), "scaffold_18", g["start"], g["end"], g["strand"], pos, g["name"], g.get("product", ""),
                    desc, pref] + flags + [note])
print(open("ahp1_region_genes.tsv").read())
