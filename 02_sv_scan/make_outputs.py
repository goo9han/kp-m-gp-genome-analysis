#!/usr/bin/env python3
"""Assemble candidates.tsv from the scan + validation JSONs, apply the verdict rules below, and draw the figure.

Interval values are means of the 1-kb bins inside the interval (the 'best window' is kept only as a column).
Verdict rules (deletion, kpgp3 coordinates):
  SUPPORTED   : not within 5 kb of a contig end; self-BLAST >= 98 % copy covers < 20 %; same-contig FR junction pairs >= 3 in each GP strain and
                0 in all 7 mouse; ONT reads skipping the interval >= 3 in Kp-GP ONT and 0 in both mouse ONT
  ARTIFACT_DUP: self-BLAST >= 98 % copy covers >= 30 % of the interval (haplotig / collapsed duplicate)
  ARTIFACT_REPEAT: mouse ONT interval depth > 3 (a multi-copy element in mouse)
  CONTIG_END  : within 5 kb of a contig end (not interpreted)
  UNRESOLVED  : depth criterion met, no junction meeting SUPPORTED
Gain: GP_ONLY_SEQUENCE (all 7 mouse interval means < 0.3: mouse reads do not map; reference-bias category, not a GP amplification);
      ARTIFACT_DUP; NOT_CONFIRMED_BY_ONT (either GP ONT interval mean < 1.3); CONTIG_END flag; GAIN (Illumina + ONT).
Reverse (kp1 coordinates, sequence missing from all GP copies): SEQUENCE_IN_GP_ASSEMBLY
      (either GP assembly covers >= 60 %: mapping, not absence); NOT_CLADE_SPECIFIC (any other mouse assembly covers < 90 %); SUPPORTED (the same junction in both GP assemblies + spanning pairs in both GP);
      else UNRESOLVED."""
import os, json, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__)); W = f"{HERE}/work"
GP = ["kpgp3_self", "kpgp6_on3"]
MOUSE = ["kp1_on3", "kpmtx_on3", "C89-2_on3", "ID1LA_on3", "IDNYC_on3", "QC1_on3", "W26G_on3"]
ONT = ["kpgp3_ont", "kpgp6_ont", "kp1_ont", "kpmtx_ont"]
NAMES = {"kpgp3_self": "Kp-GP", "kpgp6_on3": "Kp-GP2", "kp1_on3": "Kp-M", "kpmtx_on3": "Kp-MTX"}
MOUSE_C, GP_C, ACC = "#F4A261", "#3a7ebf", "#A32E2E"
HOG1_AHP1 = ("scaffold_18", 79562, 93351)

BINS = {}
for l in GP + MOUSE + ONT + ["kpgp3_on_k1", "kpgp6_on_k1", "kp1_self_k1"]:
    d = collections.defaultdict(dict)
    for r in open(f"{W}/bins_{l}.tsv"):
        c, s, e, m, n = r.rstrip("\n").split("\t"); d[c][int(s)] = float(n)
    BINS[l] = d


def imean(l, c, s, e):
    v = [x for k, x in BINS[l][c].items() if k >= s and k < e]
    return float(np.mean(v)) if v else float("nan")


def gene_str(G):
    ins = [g for g in G if g["inside"]]
    return "; ".join(f"{g['id']}{'(' + g['name'] + ')' if g['name'] else ''}" for g in ins) or "none inside"


rows = []
# ---------------- deletions ----------------
for r in json.load(open(f"{W}/validated_deletion.json")):
    c, s, e = r["contig"], r["start0"], r["end"]
    gp = [imean(g, c, s, e) for g in GP]; mo = [imean(m, c, s, e) for m in MOUSE]; on = [imean(o, c, s, e) for o in ONT]
    pr = r["pairs"]; ont = r["ont"]; sb = r["selfblast"]
    io = r.get("interchrom_ont", {})
    if r["near_contig_end"]: v = "CONTIG_END"
    elif sb["cov98"] >= 0.3: v = "ARTIFACT_DUP"
    elif max(on[2], on[3]) > 3: v = "ARTIFACT_REPEAT"
    elif (all(pr[g] >= 3 for g in GP) and all(pr[m] == 0 for m in MOUSE) and ont["kpgp3_ont"][0] >= 3
          and ont["kp1_ont"][0] == 0 and ont["kpmtx_ont"][0] == 0 and sb["cov98"] < 0.2): v = "SUPPORTED"
    else: v = "UNRESOLVED"
    note = []
    if io:
        gpt = [io[o][1] for o in ("kpgp3_ont", "kpgp6_ont")]
        if gpt[0] and gpt[0] == gpt[1] and io["kpgp3_ont"][2] >= 3 and io["kpgp6_ont"][2] >= 3 and max(io["kp1_ont"][0], io["kpmtx_ont"][0]) <= 2:
            note.append(f"GP ONT reads join this locus to {gpt[0]} (Kp-GP {io['kpgp3_ont'][2]}, Kp-GP2 {io['kpgp6_ont'][2]} reads; mouse ONT inter-contig reads "
                        f"{io['kp1_ont'][0]}/{io['kpmtx_ont'][0]})")
    if r["boundary"]["rep316_copies_near_ends"]: note.append("316-bp repeat copy within 3 kb of an end: " + ",".join(r["boundary"]["rep316_copies_near_ends"]))
    rows.append(dict(category="GP_deletion_partial_copy", coords="kpgp3", interval=f"{c}:{s + 1}-{e}", width=e - s, scales_kb=",".join(map(str, r["scales_kb"])),
                     gp_depth=f"{gp[0]:.2f}/{gp[1]:.2f}", mouse_depth=f"{min(mo):.2f}-{max(mo):.2f}", ont_gp=f"{on[0]:.2f}/{on[1]:.2f}", ont_mouse=f"{on[2]:.2f}/{on[3]:.2f}",
                     junction_pairs_gp=f"{pr['kpgp3_self']}/{pr['kpgp6_on3']}", junction_pairs_mouse=sum(pr[m] for m in MOUSE),
                     ont_skip_reads=f"GP {ont['kpgp3_ont'][0]}/{ont['kpgp6_ont'][0]}; mouse {ont['kp1_ont'][0]}/{ont['kpmtx_ont'][0]}",
                     selfblast_cov98=sb["cov98"], genes=gene_str(r["genes"]), verdict=v, note="; ".join(note),
                     is_control="YES" if c == HOG1_AHP1[0] and s < HOG1_AHP1[2] and e > HOG1_AHP1[1] else "", _gpmin=min(gp)))
order = {"SUPPORTED": 0, "UNRESOLVED": 1, "ARTIFACT_REPEAT": 2, "ARTIFACT_DUP": 3, "CONTIG_END": 4}
dels = sorted([x for x in rows], key=lambda x: (order[x["verdict"]], x["_gpmin"]))
for i, x in enumerate(dels, 1): x["rank"] = i

# ---------------- gains ----------------
gains = []
for r in json.load(open(f"{W}/validated_gain.json")):
    c, s, e = r["contig"], r["start0"], r["end"]
    gp = [imean(g, c, s, e) for g in GP]; mo = [imean(m, c, s, e) for m in MOUSE]; on = [imean(o, c, s, e) for o in ONT]
    sb = r["selfblast"]
    if max(mo) < 0.3: v = "GP_ONLY_SEQUENCE"
    elif sb["cov98"] >= 0.3: v = "ARTIFACT_DUP"
    elif min(on[0], on[1]) < 1.3: v = "NOT_CONFIRMED_BY_ONT"
    else: v = "GAIN"
    if r["near_contig_end"]: v += "|CONTIG_END"
    gains.append(dict(category="GP_gain", coords="kpgp3", interval=f"{c}:{s + 1}-{e}", width=e - s, scales_kb=",".join(map(str, r["scales_kb"])),
                      gp_depth=f"{gp[0]:.2f}/{gp[1]:.2f}", mouse_depth=f"{min(mo):.2f}-{max(mo):.2f}", ont_gp=f"{on[0]:.2f}/{on[1]:.2f}", ont_mouse=f"{on[2]:.2f}/{on[3]:.2f}",
                      junction_pairs_gp=f"RF {r['pairs']['kpgp3_self']}/{r['pairs']['kpgp6_on3']}", junction_pairs_mouse=sum(r["pairs"][m] for m in MOUSE),
                      ont_skip_reads="", selfblast_cov98=sb["cov98"], genes=gene_str(r["genes"]), verdict=v, note="", is_control="", _gpmin=min(gp)))
gains.sort(key=lambda x: (x["verdict"].split("|")[0] != "GAIN", "CONTIG_END" in x["verdict"], -x["_gpmin"]))
for i, x in enumerate(gains, 1): x["rank"] = i

# ---------------- reverse (kp1 coordinates) ----------------
junc = {j["interval"]: j for j in json.load(open(f"{W}/asm_junctions.json"))}
rev = []
for r in json.load(open(f"{W}/validated_reverse_gp_absent.json")):
    c, s, e = r["contig"], r["start0"], r["end"]; iv = f"{c}:{s + 1}-{e}"
    ac = r["assembly_coverage_ge90"]; j = junc.get(iv, {}).get("junctions", {})
    if max(ac["kpgp3"], ac["kpgp6"]) >= 0.6: v = "SEQUENCE_IN_GP_ASSEMBLY"
    elif min(ac[m] for m in ("kpmtx", "C89-2", "ID1LA", "IDNYC", "QC1", "W26G")) < 0.9: v = "NOT_CLADE_SPECIFIC"
    elif j.get("kpgp3") and j.get("kpgp6") and all(r["gp_pairs_spanning"][g] >= 3 for g in ("kpgp3_on_k1", "kpgp6_on_k1")): v = "SUPPORTED"
    else: v = "UNRESOLVED"
    pol = ""
    jj = j.get("kpgp3") or j.get("kpgp6")
    rev.append(dict(category="GP_all_copies_absent", coords="kp1", interval=iv, width=e - s, scales_kb="1",
                    gp_depth=f"{r['depth']['kpgp3_on_k1']:.2f}/{r['depth']['kpgp6_on_k1']:.2f}", mouse_depth=f"Kp-M {r['depth']['kp1_self_k1']:.2f}",
                    ont_gp="", ont_mouse="", junction_pairs_gp=f"{r['gp_pairs_spanning']['kpgp3_on_k1']}/{r['gp_pairs_spanning']['kpgp6_on_k1']}",
                    junction_pairs_mouse="", ont_skip_reads="",
                    selfblast_cov98="", genes=gene_str([dict(g, inside=True) for g in r["genes"]]), verdict=v,
                    note="; ".join(x for x in [
                        f"GP assembly coverage {ac['kpgp3']}/{ac['kpgp6']}; other mouse assemblies {min(ac[m] for m in ('kpmtx', 'C89-2', 'ID1LA', 'IDNYC', 'QC1', 'W26G')):.2f}-1.00",
                        pol,
                        (f"GP assembly junction: kp1 {jj['kp1_missing']} ({jj['missing_bp']} bp) missing; kpgp3 {j['kpgp3']['gp_contig']}:{j['kpgp3']['junction_on_gp']}"
                         + (f", kpgp6 {j['kpgp6']['gp_contig']}:{j['kpgp6']['junction_on_gp']}" if j.get("kpgp6") else "")) if jj else ""] if x),
                    is_control="", _gpmin=0))
rev.sort(key=lambda x: {"SUPPORTED": 0, "UNRESOLVED": 1, "NOT_CLADE_SPECIFIC": 2, "SEQUENCE_IN_GP_ASSEMBLY": 3}[x["verdict"]])
for i, x in enumerate(rev, 1): x["rank"] = i

# ---------------- GP-present / mouse-absent (reference-bias category) ----------------
mu = []
for r in json.load(open(f"{W}/validated_mouse_unmappable.json")):
    c, s, e = r["contig"], r["start0"], r["end"]; d = r["depth"]
    sk = r["mouse_ont_skip"]
    mu.append(dict(category="GP_present_mouse_absent", coords="kpgp3", interval=f"{c}:{s + 1}-{e}", width=e - s, scales_kb="1",
                   gp_depth=f"{d['kpgp3_self']:.2f}/{d['kpgp6_on3']:.2f}", mouse_depth=f"{min(d[m] for m in MOUSE):.2f}-{max(d[m] for m in MOUSE):.2f}",
                   ont_gp=f"{d['kpgp3_ont']:.2f}/{d['kpgp6_ont']:.2f}", ont_mouse=f"{d['kp1_ont']:.2f}/{d['kpmtx_ont']:.2f}", junction_pairs_gp="",
                   junction_pairs_mouse="", ont_skip_reads=f"mouse ONT reads spanning the locus without it: {sk['kp1_ont']}/{sk['kpmtx_ont']}",
                   selfblast_cov98="", genes=gene_str(r["genes"]),
                   verdict="MOUSE_READS_SPAN_LOCUS" if min(sk.values()) >= 3 else ("CONTIG_END" if r["near_contig_end"] else "NO_SPANNING_MOUSE_READS"),
                   note="", is_control="", _gpmin=0))
mu.sort(key=lambda x: (x["verdict"] != "MOUSE_READS_SPAN_LOCUS", -x["width"]))
for i, x in enumerate(mu, 1): x["rank"] = i

cols = ["category", "rank", "coords", "interval", "width", "scales_kb", "gp_depth", "mouse_depth", "ont_gp", "ont_mouse", "junction_pairs_gp",
        "junction_pairs_mouse", "ont_skip_reads", "selfblast_cov98", "genes", "verdict", "note", "is_control"]
with open(f"{HERE}/candidates.tsv", "w") as fh:
    fh.write("# sv_scan candidates. Depth = mean of 1-kb bins inside the interval, each strain normalized to its own genome median (MAPQ >= 20).\n"
             "# gp_depth = Kp-GP/Kp-GP2; mouse_depth = min-max of the 7 mouse strains; ont = Kp-GP/Kp-GP2 and Kp-M/Kp-MTX ONT.\n"
             "# junction_pairs = unique FR (deletion) or RF (gain) Illumina pairs spanning the interval; categories and verdict rules: make_outputs.py docstring.\n"
             "# NOTE: the 'interval' column is the SCAN WINDOW the scan called (window-level edges), not the validated deletion boundary.\n"
             "#       Validated boundaries: HOG1-AHP1 scaffold_18:79562-93351 (13.8 kb); CHK1 scaffold_6:420988-422860 (1.9 kb). See figures/panels/svscan_numbers.md.\n")
    fh.write("\t".join(cols) + "\n")
    for x in dels + gains + rev + mu: fh.write("\t".join(str(x[k]) for k in cols) + "\n")
json.dump(dict(deletion=dels, gain=gains, reverse=rev, mouse_absent=mu), open(f"{W}/final.json", "w"), indent=1, default=str)
print("deletions:", collections.Counter(x["verdict"] for x in dels)); print("gains:", collections.Counter(x["verdict"] for x in gains))
print("reverse:", collections.Counter(x["verdict"] for x in rev)); print("GP-present/mouse-absent:", collections.Counter(x["verdict"] for x in mu), sum(x["width"] for x in mu), "bp")
print("HOG1-AHP1 rank:", [x["rank"] for x in dels if x["is_control"]])

# ---------------- figure ----------------
plt.rcParams.update({"font.family": "Arial", "svg.fonttype": "none", "font.size": 7, "text.color": "#222", "axes.edgecolor": "#666",
                     "xtick.color": "#666", "ytick.color": "#666"})
k3 = BINS["kpgp3_self"]
ctgs = sorted(k3, key=lambda c: -max(k3[c]))
off, x0 = {}, 0
for c in ctgs: off[c] = x0; x0 += max(k3[c]) + 1000 + 30000
fig = plt.figure(figsize=(7.2, 5.2))
gs = fig.add_gridspec(3, 4, height_ratios=[1.0, 1.25, 1.35], hspace=0.55, wspace=0.35)
ax1 = fig.add_subplot(gs[0, :]); ax2 = fig.add_subplot(gs[1, :], sharex=ax1)
X, gpv, mov = [], [], []
for c in ctgs:
    ks = sorted(k3[c])
    for i in range(0, len(ks) - 4, 5):   # tiled 5-kb windows
        b = ks[i:i + 5]
        X.append(off[c] + b[0] + 2500)
        gpv.append(np.mean([np.mean([BINS[g][c][k] for k in b]) for g in GP]))
        mov.append(np.mean([np.mean([BINS[m][c][k] for k in b]) for m in MOUSE]))
X, gpv, mov = np.array(X), np.array(gpv), np.array(mov)
for i, c in enumerate(ctgs):
    if i % 2 == 0:
        for ax in (ax1, ax2): ax.axvspan(off[c] / 1e6, (off[c] + max(k3[c]) + 1000) / 1e6, color="#f3f4f6", lw=0, zorder=0)
ax1.scatter(X / 1e6, np.clip(mov, 0, 3), s=0.6, color=MOUSE_C, lw=0, label="mouse (mean of 7)", rasterized=True)
ax1.scatter(X / 1e6, np.clip(gpv, 0, 3), s=0.6, color=GP_C, lw=0, label="guinea pig (mean of Kp-GP, Kp-GP2)", rasterized=True)
ax1.set_ylim(0, 3.05); ax1.set_yticks([0, 1, 2, 3]); ax1.set_ylabel("normalized depth\n(5-kb windows)")
ax1.legend(loc="upper right", fontsize=6, frameon=False, markerscale=6, ncol=2, handletextpad=0.2)
ax1.text(0, 3.3, "A", fontsize=10, fontweight="bold", va="bottom")
ok = mov >= 0.3
ratio = np.log2(np.clip(gpv, 0.02, None) / np.clip(mov, 0.02, None))
ax2.scatter(X[ok] / 1e6, np.clip(ratio[ok], -2.5, 2.5), s=0.6, color="#6b7280", lw=0, rasterized=True)
ax2.scatter(X[~ok] / 1e6, np.full((~ok).sum(), 2.45), s=2, marker="|", color="#9ca3af", lw=0.5, rasterized=True)
ax2.axhline(0, color="#999", lw=0.5)
for y in (np.log2(2 / 3), np.log2(1 / 3)): ax2.axhline(y, color="#ccc", lw=0.4, ls=(0, (2, 2)))
ax2.set_ylim(-2.7, 2.7); ax2.set_ylabel("log2(GP / mouse)")
ax2.text(X.max() / 1e6, 2.62, "ticks at top: mouse < 0.3 (GP-only sequence, ratio not shown)", ha="right", va="bottom", fontsize=5.5, color="#6b7280")
mk = {"SUPPORTED": (ACC, "v", 22), "UNRESOLVED": ("#b45309", "v", 12)}
for x in dels:
    if x["verdict"] in mk:
        c, a = x["interval"].split(":"); s, e = map(int, a.split("-"))
        col, m, sz = mk[x["verdict"]]
        ax2.scatter([(off[c] + (s + e) / 2) / 1e6], [-2.45], marker="^", s=sz, color=col, zorder=5)
for x in gains:
    if x["verdict"].split("|")[0] == "GAIN":
        c, a = x["interval"].split(":"); s, e = map(int, a.split("-"))
        ax2.scatter([(off[c] + (s + e) / 2) / 1e6], [2.2], marker="v", s=14, color=GP_C, zorder=5)
ax2.set_xlabel("Kp-GP assembly (contigs >= 20 kb, longest first; Mb)")
ax2.text(0, 2.95, "B", fontsize=10, fontweight="bold", va="bottom")
ax2.scatter([], [], marker="^", s=22, color=ACC, label="GP deletion, junction-supported")
ax2.scatter([], [], marker="^", s=12, color="#b45309", label="GP depth drop, junction unresolved")
ax2.scatter([], [], marker="v", s=14, color=GP_C, label="GP gain (Illumina + ONT)")
ax2.legend(loc="lower right", fontsize=5.8, frameon=False, ncol=3, handletextpad=0.2, bbox_to_anchor=(1.0, 0.02))
# zooms
Z = [("scaffold_18", 70000, 105000, "HOG1–AHP1 (Kp-GP scaffold_18)", "k3", (79562, 93351)),
     ("scaffold_3", 85000, 115000, "Kp-GP scaffold_3:97–104 kb", "k3", (97000, 104000)),
     ("scaffold_55", 0, 47936, "Kp-GP scaffold_55 (whole contig)", "k3", None),
     ("scaffold_6", 405000, 440000, "CHK1 region (Kp-M scaffold_6)", "k1", (420988, 422860))]
for i, (c, a, b, t, tag, hl) in enumerate(Z):
    ax = fig.add_subplot(gs[2, i])
    if tag == "k3":
        for g, ls in zip(GP, ("-", (0, (2, 1)))):
            ks = sorted(k for k in BINS[g][c] if a <= k < b); ax.plot([(k + 500) / 1e3 for k in ks], [BINS[g][c][k] for k in ks], color=GP_C, lw=0.8, ls=ls)
        for m in MOUSE:
            ks = sorted(k for k in BINS[m][c] if a <= k < b); ax.plot([(k + 500) / 1e3 for k in ks], [BINS[m][c][k] for k in ks], color=MOUSE_C, lw=0.5, alpha=0.7)
    else:
        for g, ls in zip(("kpgp3_on_k1", "kpgp6_on_k1"), ("-", (0, (2, 1)))):
            ks = sorted(k for k in BINS[g][c] if a <= k < b); ax.plot([(k + 500) / 1e3 for k in ks], [BINS[g][c][k] for k in ks], color=GP_C, lw=0.8, ls=ls)
        ks = sorted(k for k in BINS["kp1_self_k1"][c] if a <= k < b)
        ax.plot([(k + 500) / 1e3 for k in ks], [BINS["kp1_self_k1"][c][k] for k in ks], color=MOUSE_C, lw=0.8)
    if hl: ax.axvspan(hl[0] / 1e3, hl[1] / 1e3, color="#FDECEC", lw=0, zorder=0)
    ax.axhline(1, color="#bbb", lw=0.4)
    ax.set_ylim(0, 3.0 if c == "scaffold_55" else 1.8); ax.set_title(t, fontsize=6.3, pad=3)
    ax.set_xlabel("kb", fontsize=6); ax.tick_params(labelsize=5.5)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    if i == 0:
        ax.set_ylabel("normalized depth", fontsize=6.3); ax.text(a / 1e3, 2.05, "C", fontsize=10, fontweight="bold", va="bottom")
for ax in (ax1, ax2):
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
plt.setp(ax1.get_xticklabels(), visible=False)
fig.subplots_adjust(left=0.08, right=0.99, top=0.95, bottom=0.08)
for ext in ("svg", "png"): fig.savefig(f"{HERE}/sv_scan_genome.{ext}", dpi=300)
print("figure saved")
