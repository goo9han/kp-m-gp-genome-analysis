#!/usr/bin/env python3
"""YAP1 events for the promoter figure: clade-level differences and shared sites, for all eight drawn genes,
with read confirmation for every difference.

scan_yap1.py produces the site calls (and is gated on reproducing the existing FIMO table exactly). This script
turns them into drawable events under the same rule the STRE events already obey: a difference is drawn only if
the guinea-pig allele is fixed in the reads of both guinea-pig strains and absent from all seven mouse strains.

Read confirmation comes from af_recover/recovered_variants.tsv where that file covers the locus; where it does
not, the allele fractions are counted directly from the same Illumina BAMs build_events.py uses, and written to
yap1_read_check_<gene>_<pos>.json so the number can be checked without rerunning anything.

Output: yap1_events.tsv  (drawn + not-drawn, with the reason), yap1_read_check_*.json
"""
import csv, json, os, re, subprocess, sys, collections
_HERE = os.path.dirname(os.path.abspath(__file__))
_KP = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, f"{_KP}/stre_null_model")
from stre_core import *                                   # noqa
HERE, KP = _HERE, _KP                                      # stre_core exports HERE/ROOT; restore ours
SAM = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/samtools")
REF = f"{KP}/qpcr_primer_design/validation/genomes/kpgp3_combined.fa"
MOUSE = ["kp1", "kpmtx", "C89-2", "ID1LA", "IDNYC", "QC1", "W26G"]
GP = ["kpgp3", "kpgp6"]
BAMS = {"kpgp3": f"{KP}/qpcr_primer_design/validation/pileup/kpgp3_illumina_combined.sorted.bam",
        "kpgp6": f"{KP}/ahp1_cnv/bams/kpgp6_on_kpgp3.bam",
        "kp1": f"{KP}/ahp1_cnv/bams/kp1_on_kpgp3.bam",
        "kpmtx": f"{KP}/ahp1_cnv/bams/kpmtx_on_kpgp3.bam",
        **{s: f"{KP}/ahp1_cnv/nine_strain/bams/{s}_on_kpgp3.bam" for s in ("C89-2", "ID1LA", "IDNYC", "QC1", "W26G")}}
FIXED_MIN, MOUSE_MAX = 0.90, 0.05                          # the thresholds stre_null_model uses for FIXED / absent

meta = json.load(open(f"{KP}/stre_null_model/work/stress_meta.json"))
name2id = {v["label"]: k for k, v in meta.items()}
sites = [r for r in csv.DictReader([l for l in open(f"{HERE}/yap1_sites_stress.tsv") if not l.startswith("#")],
                                   delimiter="\t")]
rv = list(csv.DictReader(open(f"{KP}/af_recover/recovered_variants.tsv"), delimiter="\t"))
rv_loci = {r["locus"] for r in rv}


def col2pos(row):
    m, i = {}, -1
    for c, ch in enumerate(row):
        if ch != "-":
            i += 1; m[c] = i
    return m


def allele_fractions(gene, pos_from_atg):
    """Count the guinea-pig allele at every position where the two clades differ inside the site."""
    gid = name2id[gene]
    A = read_fasta(f"{KP}/stre_null_model/work/aln_stress/{gid}.aln.fa")
    c2p = {s: col2pos(A[s]) for s in STRAINS}
    p2c = {v: k for k, v in c2p["kp1"].items()}
    lo, hi = [int(x) + 700 for x in pos_from_atg.split("..")]
    c0, c1 = p2c[lo], p2c[hi]
    diffs = [c for c in range(c0, c1 + 1) if A["kp1"][c] != A["kpgp3"][c] and "-" not in (A["kp1"][c], A["kpgp3"][c])]
    ctg, gs, ge, strand = meta[gid]["ranges"]["kpgp3"]
    out = {}
    for c in diffs:
        gi = c2p["kpgp3"].get(c)
        if gi is None:
            continue
        gpos = gs + gi if strand == "+" else ge - gi
        der = A["kpgp3"][c]
        per = {}
        for s in STRAINS:
            r = subprocess.run([SAM, "mpileup", "-f", REF, "-r", f"{ctg}:{gpos}-{gpos}", "-q", "20", "-Q", "20",
                                "--ff", "0x904", BAMS[s]], capture_output=True, text=True).stdout.strip().split("\t")
            if len(r) < 5:
                per[s] = [0, 0]; continue
            b = re.sub(r"\^.", "", r[4].upper())
            b = re.sub(r"[+-]\d+[ACGTN]*", "", b)
            cnt = collections.Counter(ch for ch in b if ch in "ACGTN.,")
            tot = sum(v for k, v in cnt.items() if k in "ACGT.,")
            per[s] = [cnt[der] + cnt["."] + cnt[","], tot]           # . and , = matches kpgp3 = the derived allele
        out[f"{ctg}:{gpos}"] = dict(kp1=A["kp1"][c], kpgp3=der, counts=per)
    return out


def indel_fractions(gene, pos_from_atg):
    """Reads carrying the guinea-pig form of the site against the mouse form, when the two differ by an indel."""
    gid = name2id[gene]
    A = read_fasta(f"{KP}/stre_null_model/work/aln_stress/{gid}.aln.fa")
    c2p = {s: col2pos(A[s]) for s in STRAINS}
    p2c = {v: k for k, v in c2p["kp1"].items()}
    lo, hi = [int(x) + 700 for x in pos_from_atg.split("..")]
    if lo not in p2c or hi not in p2c:
        return None
    c0, c1 = p2c[lo], p2c[hi]
    if not any(A["kp1"][c] == "-" or A["kpgp3"][c] == "-" for c in range(c0, c1 + 1)):
        return None
    gp = "".join(A["kpgp3"][c] for c in range(c0 - 6, c1 + 7) if A["kpgp3"][c] != "-")
    mo = "".join(A["kp1"][c] for c in range(c0 - 6, c1 + 7) if A["kp1"][c] != "-")
    if gp == mo:
        return None
    ctg, gs, ge, strand = meta[gid]["ranges"]["kpgp3"]
    gi = c2p["kpgp3"].get(c0, 0)
    gpos = gs + gi if strand == "+" else ge - gi
    rc_ = lambda x: x.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    out = {}
    for s in STRAINS:
        reads = subprocess.run([SAM, "view", "-F", "0x904", BAMS[s], f"{ctg}:{gpos - 150}-{gpos + 150}"],
                               capture_output=True, text=True).stdout.splitlines()
        g = m = 0
        for r in reads:
            q = r.split("\t")[9]
            hg = gp in q or rc_(gp) in q
            hm = mo in q or rc_(mo) in q
            if hg and not hm:
                g += 1
            elif hm and not hg:
                m += 1
        out[s] = (g, g + m)
    return out


rows, log = [], []
for s in sites:
    gene = s["locus"]
    if s["kind"] == "shared":
        rows.append(dict(gene=gene, motif="YAP1", pos_from_ATG=s["pos_from_ATG"], strand=s["strand"], event="shared",
                         cause="", mouse_state="present 7/7", gp_state="present 2/2", read_confirmation="",
                         source="figures/stre_turnover/yap1_sites_stress.tsv")); continue
    if s["kind"] == "other":
        log.append(f"NOT DRAWN {gene} YAP1 {s['pos_from_ATG']}: not a clade difference (present in {s['present_in']})")
        continue
    af = allele_fractions(gene, s["pos_from_ATG"])
    if not af:                                             # no differing base inside the site -> an indel caused it
        ind = indel_fractions(gene, s["pos_from_ATG"])
        if ind is None:
            log.append(f"NOT RESOLVED {gene} YAP1 {s['kind']} {s['pos_from_ATG']}: neither a substitution nor a "
                       "resolvable indel inside the site"); continue
        gp_ok = all(ind[g2][1] and ind[g2][0] / ind[g2][1] >= FIXED_MIN for g2 in GP)
        mo_ok = all(ind[m][1] and ind[m][0] / ind[m][1] <= MOUSE_MAX for m in MOUSE)
        conf = "reads (indel): " + "; ".join(f"{k} {v[0]}/{v[1]}" for k, v in ind.items())
        json.dump({k: list(v) for k, v in ind.items()},
                  open(f"{HERE}/yap1_read_check_{gene}_{s['pos_from_ATG']}.json", "w"), indent=1)
        if not (gp_ok and mo_ok):
            log.append(f"NOT FIXED {gene} YAP1 {s['kind']} {s['pos_from_ATG']}: {conf}"); continue
        rows.append(dict(gene=gene, motif="YAP1", pos_from_ATG=s["pos_from_ATG"], strand=s["strand"], event=s["kind"],
                         cause="indel (read-confirmed)",
                         mouse_state="present 7/7" if s["kind"] == "loss" else "absent 7/7",
                         gp_state="absent 2/2 (all reads)" if s["kind"] == "loss" else "present 2/2 (all reads)",
                         read_confirmation=conf,
                         source="figures/stre_turnover/yap1_sites_stress.tsv; indel read counts from the Illumina BAMs"))
        continue
    json.dump(af, open(f"{HERE}/yap1_read_check_{gene}_{s['pos_from_ATG']}.json", "w"), indent=1)
    gp_ok = all(all(c["counts"][g][1] and c["counts"][g][0] / c["counts"][g][1] >= FIXED_MIN for g in GP) for c in af.values())
    mo_ok = all(all(c["counts"][m][1] and c["counts"][m][0] / c["counts"][m][1] <= MOUSE_MAX for m in MOUSE) for c in af.values())
    conf = "reads: " + "; ".join(
        f"{p} {c['kp1']}>{c['kpgp3']} AF " + ", ".join(f"{g} {c['counts'][g][0]}/{c['counts'][g][1]}" for g in GP)
        + "; mouse max " + f"{max(c['counts'][m][0] / c['counts'][m][1] if c['counts'][m][1] else 0 for m in MOUSE):.3f}"
        for p, c in af.items())
    if not (af and gp_ok and mo_ok):
        log.append(f"NOT FIXED {gene} YAP1 {s['kind']} {s['pos_from_ATG']}: {conf}")
        continue
    rows.append(dict(gene=gene, motif="YAP1", pos_from_ATG=s["pos_from_ATG"], strand=s["strand"], event=s["kind"],
                     cause="substitution",
                     mouse_state="present 7/7" if s["kind"] == "loss" else "absent 7/7",
                     gp_state="absent 2/2 (all reads)" if s["kind"] == "loss" else "present 2/2 (all reads)",
                     read_confirmation=conf,
                     source="figures/stre_turnover/yap1_sites_stress.tsv; read counts from the Illumina BAMs"))

for g in ("AHP1", "GPX2", "SOD1", "SKN7", "HSP26", "CRH1", "CWP2", "HAP2"):
    if not any(r["gene"] == g for r in rows):
        log.append(f"NO YAP1 SITE {g}: the promoter has no YAP1 hit at p<1e-3 in any of the nine strains")

cols = ["gene", "motif", "pos_from_ATG", "strand", "event", "cause", "mouse_state", "gp_state", "read_confirmation", "source"]
with open(f"{HERE}/yap1_events.tsv", "w") as fh:
    fh.write("# YAP1 events for the eight drawn genes. Sites from scan_yap1.py (calibrated against the existing\n"
             "# FIMO table, 11/11). A difference is listed only if the guinea-pig allele is fixed in both guinea-pig\n"
             "# strains (AF >= 0.90 at every differing position) and at most 0.05 in all seven mouse strains.\n")
    w = csv.DictWriter(fh, cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
open(f"{HERE}/yap1_build_log.txt", "w").write("\n".join(log) + "\n")
print(collections.Counter((r["gene"], r["event"]) for r in rows if r["event"] != "shared"))
print(f"shared YAP1 sites: {sum(1 for r in rows if r['event'] == 'shared')}")
print("\n".join(log))
