#!/usr/bin/env python3
"""Event table for the promoter cis-element figure (STRE + YAP1, GP-specific FIXED changes).

No motif scanning or AF calling is redone here. Sources:
  STRE events      stre_null_model/stre_events_filtered.tsv   (status PASS = all causal substitutions FIXED in both GP strains;
                                                                EXCLUDED_FROM_TESTS = indel, kept only if the read check below confirms it)
  STRE context     stre_turnover_figure/stre_turnover_data.tsv (sites shared by all 9 strains; same site calls as the events file)
  YAP1 sites       h2o2_analysis/true_atg_all/new/sites_all.tsv (true-ATG promoters of the 14 H2O2 loci only)
  YAP1 read AF     af_recover/recovered_variants.tsv            (causal substitution with 'MA0415.2 gain/loss', joint_class FIXED)
The only new computation is a read-level homopolymer length count for indel events (samtools view on existing BAMs).
Positions are relative to the kp1 ATG (alignment-anchored, -700..-1).
"""
import csv, re, subprocess, collections, json, os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.dirname(os.path.abspath(__file__))
SAMTOOLS = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/samtools")
REF = f"{ROOT}/qpcr_primer_design/validation/genomes/kpgp3_combined.fa"
MOUSE = ["kp1", "kpmtx", "C89-2", "ID1LA", "IDNYC", "QC1", "W26G"]
GP = ["kpgp3", "kpgp6"]
# Illumina BAMs on the kpgp3 assembly (existing; see ahp1_cnv/ and nine_strain/)
BAMS = {"kpgp3": f"{ROOT}/qpcr_primer_design/validation/pileup/kpgp3_illumina_combined.sorted.bam",
        "kpgp6": f"{ROOT}/ahp1_cnv/bams/kpgp6_on_kpgp3.bam",
        "kp1": f"{ROOT}/ahp1_cnv/bams/kp1_on_kpgp3.bam",
        "kpmtx": f"{ROOT}/ahp1_cnv/bams/kpmtx_on_kpgp3.bam",
        **{s: f"{ROOT}/ahp1_cnv/nine_strain/bams/{s}_on_kpgp3.bam" for s in ("C89-2", "ID1LA", "IDNYC", "QC1", "W26G")}}

# Genes named in the brief, plus HAP2 (meets the same criterion: FIXED STRE event in stre_events_filtered.tsv; a gain).
REQUESTED = ["AHP1", "SOD1", "SKN7", "CTA1", "GPX2", "TRX2", "HSP26", "CRH1", "CWP2"]
EXTRA = ["HAP2"]
OXIDATIVE = ["AHP1", "GPX2", "SOD1", "SKN7", "CTA1", "TRX2"]   # the loci the old H2O2-set YAP1 scan covered; kept for the log only
log = []


def read_tsv(p, skip_comments=True):
    with open(p) as fh:
        rows = [l for l in fh if not (skip_comments and l.startswith("#"))]
    return list(csv.DictReader(rows, delimiter="\t"))


def homopolymer_check(contig, pos):
    """Count homopolymer run lengths in reads of all 9 strains at an indel site (kpgp3 coordinates)."""
    seq = subprocess.run([SAMTOOLS, "faidx", REF, f"{contig}:{pos - 20}-{pos + 20}"], capture_output=True, text=True).stdout
    seq = "".join(seq.split("\n")[1:]).upper()
    runs = [(m.start(), m.end()) for m in re.finditer(r"(A+|C+|G+|T+)", seq) if m.end() - m.start() >= 3]
    a, b = min(runs, key=lambda r: abs((r[0] + r[1]) / 2 - 20))
    base, left, right = seq[a], seq[a - 8:a], seq[b:b + 8]
    pat = re.compile(left + f"({base}+)" + right)
    res = {}
    for s, bam in BAMS.items():
        reads = subprocess.run([SAMTOOLS, "view", "-F", "0x904", bam, f"{contig}:{pos - 5}-{pos + 5}"],
                               capture_output=True, text=True).stdout.splitlines()
        c = collections.Counter()
        for r in reads:
            m = pat.search(r.split("\t")[9])
            c[len(m.group(1)) if m else "no_flank_match"] += 1
        res[s] = dict(c)
    return dict(ref_context=f"{left}[{base}x{b - a}]{right}", counts=res)


ev = [r for r in read_tsv(f"{ROOT}/stre_null_model/stre_events_filtered.tsv") if r["set"] == "STRESS"]
ctx = read_tsv(f"{ROOT}/stre_turnover_figure/stre_turnover_data.tsv")
sites = read_tsv(f"{ROOT}/h2o2_analysis/true_atg_all/new/sites_all.tsv", skip_comments=False)
rv = read_tsv(f"{ROOT}/af_recover/recovered_variants.tsv", skip_comments=False)
indel_json = json.load(open(f"{ROOT}/stre_null_model/work/events_STRESS.json"))

rows = []
# ---- STRE events ----
for r in ev:
    g = r["locus"]
    kind = r["event"]
    if r["status"] == "PASS":
        conf = "; ".join(re.findall(r"[ACGT]>[ACGT]@-?\d+ \[FIXED; AF kpgp3 [\d.]+\(d=\d+\), kpgp6 [\d.]+\(d=\d+\)", r["causal_substitutions"]))
        conf = "reads: " + conf.replace(" [FIXED; ", " ") + " (derived-allele fraction in GP reads)"
        cause = "substitution"
    elif r["status"] == "EXCLUDED_FROM_TESTS" and r["cause"] == "indel":
        o = [o for o in indel_json[r["kp1_gene"]]["obs"] if o["kp1_pos"] == r["pos_from_ATG(kp1)"]][0]
        contig, pos, _ = o["gcoord"]["kpgp3"]
        hc = homopolymer_check(contig, pos)
        gp_len = {collections.Counter({k: v for k, v in hc["counts"][s].items() if k != "no_flank_match"}).most_common(1)[0][0] for s in GP}
        mo_len = {collections.Counter({k: v for k, v in hc["counts"][s].items() if k != "no_flank_match"}).most_common(1)[0][0] for s in MOUSE}
        gp_pure = all(len([k for k in hc["counts"][s] if k != "no_flank_match"]) == 1 for s in GP + MOUSE)
        ok = len(gp_len) == 1 and len(mo_len) == 1 and gp_len != mo_len and gp_pure
        json.dump(hc, open(f"{OUT}/indel_read_check_{r['locus']}_{r['pos_from_ATG(kp1)']}.json", "w"), indent=1)
        conf = ("reads: homopolymer " + hc["ref_context"] + " (kpgp3 " + f"{contig}:{pos}" + "); run length in every read "
                + ", ".join(f"{s} {list(hc['counts'][s].items())}" for s in GP + MOUSE))
        if not ok:
            log.append(f"DROPPED {g} {r['pos_from_ATG(kp1)']} indel: read check not clean ({conf})")
            continue
        cause = "1-bp indel (homopolymer)"
    else:
        log.append(f"NOT FIXED {g} STRE {kind} {r['pos_from_ATG(kp1)']}: {r['status']} {r['causal_substitutions']}")
        continue
    mouse_state = "present 7/7" if kind == "loss" else "absent 7/7"
    gp_state = "absent 2/2 (all reads)" if kind == "loss" else "present 2/2 (all reads)"
    rows.append(dict(gene=g, motif="STRE", pos_from_ATG=r["pos_from_ATG(kp1)"], strand=r["strand"], event=kind,
                     mouse_state=mouse_state, gp_state=gp_state, seq_mouse=r["seven_strain_site_seq"].split("x")[0],
                     seq_gp=r["kpgp3_site_seq"], cause=cause, read_confirmation=conf,
                     source="stre_null_model/stre_events_filtered.tsv"))

# ---- YAP1 events (all eight genes; scan_yap1.py -> yap1_events.py) ----
# yap1_sites_stress.tsv covers every gene in the stress set, and every clade difference below carries its
# own read check.
for r in read_tsv(f"{OUT}/yap1_events.tsv"):
    if r["event"] == "shared":
        continue
    rows.append(dict(gene=r["gene"], motif="YAP1", pos_from_ATG=r["pos_from_ATG"], strand=r["strand"],
                     event=r["event"], mouse_state=r["mouse_state"], gp_state=r["gp_state"],
                     seq_mouse="", seq_gp="", cause=r["cause"], read_confirmation=r["read_confirmation"],
                     source=r["source"]))
# YAP1 changes in the requested genes that are not fixed (reported, not drawn)
for v in rv:
    if v["locus"] in REQUESTED and v["where"].startswith("promoter ") and "MA0415.2" in v["motif_effect"] and v["joint_class"] != "FIXED":
        log.append(f"NOT FIXED {v['locus']} YAP1 {v['motif_effect']} {v['where']} {v['change']}: {v['joint_class']} "
                   f"AF kpgp3 {v['kpgp3_AF(der/depth)']}, kpgp6 {v['kpgp6_AF(der/depth)']}")

drawn = sorted({r["gene"] for r in rows})
# ---- context: sites shared by all 9 strains ----
for c in ctx:
    if c["locus"] in drawn and c["class"] == "shared":
        rows.append(dict(gene=c["locus"], motif="STRE", pos_from_ATG=c["pos"], strand=c["strand"], event="shared",
                         mouse_state="present 7/7", gp_state="present 2/2", seq_mouse=c["seq_kp1"], seq_gp=c["seq_kpgp3"],
                         cause="", read_confirmation="", source="stre_turnover_figure/stre_turnover_data.tsv"))
    elif c["locus"] in drawn and c["class"] == "other":
        log.append(f"NOT DRAWN {c['locus']} STRE {c['pos']}: polymorphic within a clade (present in {c['present_in']})")
for r in read_tsv(f"{OUT}/yap1_events.tsv"):
    if r["event"] == "shared" and r["gene"] in drawn:
        rows.append(dict(gene=r["gene"], motif="YAP1", pos_from_ATG=r["pos_from_ATG"], strand=r["strand"], event="shared",
                         mouse_state="present 7/7", gp_state="present 2/2", seq_mouse="", seq_gp="",
                         cause="", read_confirmation="", source=r["source"]))

for g in REQUESTED:
    if g not in drawn:
        log.append(f"GENE NOT DRAWN {g}: no GP-specific FIXED STRE/YAP1 change in the sources")
for g in EXTRA:
    if g in drawn:
        log.append(f"GENE ADDED {g}: not in the brief's list but has a FIXED STRE event in stre_events_filtered.tsv (a gain); "
                   "drawn so that gains and losses are shown under the same rule")
# fixed motif changes outside the site-level scans (reported, not drawn)
for v in rv:
    if (v["where"].startswith("promoter ") and v["joint_class"] == "FIXED" and re.search(r"MA0415|MA0341|MA0342", v["motif_effect"])
            and v["locus"] not in drawn):
        log.append(f"OUTSIDE SITE-LEVEL SCAN {v['locus']} {v['where']} {v['change']} {v['motif_effect']} (single-variant PWM effect, "
                   f"af_recover; AF kpgp3 {v['kpgp3_AF(der/depth)']}, kpgp6 {v['kpgp6_AF(der/depth)']}); not drawn")

for r in rows:
    r["yap1_scanned"] = "yes"          # all eight, since scan_yap1.py
    r["within_HOG1_AHP1_deletion"] = "yes (site change is in the retained copies)" if r["gene"] == "AHP1" else "no"
order = {"loss": 0, "gain": 1, "shared": 2}
rows.sort(key=lambda r: (r["gene"], order[r["event"]], int(r["pos_from_ATG"].split("..")[0])))
cols = ["gene", "motif", "pos_from_ATG", "strand", "event", "mouse_state", "gp_state", "seq_mouse", "seq_gp", "cause",
        "read_confirmation", "yap1_scanned", "within_HOG1_AHP1_deletion", "source"]
with open(f"{OUT}/stre_turnover_events.tsv", "w") as fh:
    fh.write("# GP-specific promoter STRE/YAP1 changes (FIXED in both GP strains, absent/present in all 7 mouse strains) and shared context sites.\n"
             "# Positions relative to the kp1 ATG (alignment-anchored). event: loss / gain / shared (all 9 strains). Built by build_events.py.\n")
    w = csv.DictWriter(fh, cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
with open(f"{OUT}/build_log.txt", "w") as fh:
    fh.write("\n".join(log) + "\n")
print(collections.Counter((r["gene"], r["motif"], r["event"]) for r in rows if r["event"] != "shared"))
print("\n".join(log))
