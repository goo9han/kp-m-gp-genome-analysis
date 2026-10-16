#!/usr/bin/env python3
"""hog1_ahp1_deletion_figure.svg/png - HOG1-AHP1 13.8-kb deletion in the guinea-pig clade.
Only existing results are read (no re-estimation):
  depth      ahp1_cnv/profile_{kpgp3_self,kpgp6_on3,kp1_on3,kpmtx_on3,kpgp3_ont,kpgp6_ont}.json (100-bp bins, MAPQ>=20,
             each strain normalised to its own genome-wide median depth on the kpgp3 assembly)
  genes      ahp1_cnv/ahp1_region_genes.tsv (kpgp3 funannotate)
  pairs      ahp1_cnv/evidence_{kpgp3_self,kpgp6_on3,kp1_on3,kpmtx_on3}.json (discordant pairs)
  ONT reads  ont_reads_{kpgp3,kpgp6}.json (extents of the deletion reads of ahp1_cnv/ont_*.json; verified identical)
Style: Arial, text #222 / axes #666, deletion accent #A32E2E on #FDECEC (as h2o2_analysis/figure_sequences figures);
clade colours mouse #1F5C99, guinea pig #D9730D. SVG text kept as text (svg.fonttype none)."""
import json, csv, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, Arc
plt.rcParams.update({"svg.fonttype": "none", "font.family": "Arial", "font.size": 7.5, "axes.edgecolor": "#666",
                     "xtick.color": "#666", "ytick.color": "#666", "axes.labelcolor": "#222", "text.color": "#222"})
HERE = os.path.dirname(os.path.abspath(__file__)); CNV = os.path.dirname(HERE)

# ---------------- naming (single source for every label) ----------------
STRAIN_NAMES = {"kp1": "Kp-M", "kpmtx": "Kp-MTX", "kpgp3": "Kp-GP", "kpgp6": "Kp-GP2"}      # fixed by the user
GENE_NAMES = {  # kpgp3 funannotate id -> label (Sc ortholog names from eggNOG / funannotate; edit here only)
    "FUN_002902": "AVL9", "FUN_002903": "tRNA-Arg", "FUN_002904": "HOG1", "FUN_002905": "CCW12", "FUN_002906": "AHP1",
    "FUN_002907": "FUN_002907", "FUN_002908": "tRNA-Trp"}
BOLD = {"FUN_002904", "FUN_002906"}
MOUSE_C, GP_C, ACC, ACC_F = "#F4A261", "#3a7ebf", "#A32E2E", "#FDECEC"   # mouse orange, guinea pig blue (as the HSP26/AHP1 promoter figures)
X0, X1 = 76000, 97000
DEL = (79562, 93351)                     # deleted interval (one repeat copy retained)
REP = [(79562, 79877), (93352, 93667)]   # the two 316-bp copies of the flanking repeat
MOUSE_ABSENT = (80500, 82300)            # blastn: kpgp3 80,500-82,300 not present at the kp1 locus (ahp1_cnv_verdict.md section 6)
MOUSE_ZERO = (79601, 82400)              # both mouse strains' plotted depth <= 0.27 (read from the plotted bins, section 1 of
                                         # deletion_figure_data.tsv); flanking bins 79,501 (0.72) and 82,401 (0.58 / 0.86)
BAND = "#e5e7eb"                         # neutral grey, distinct from the deletion accent

prof = {k: json.load(open(f"{CNV}/profile_{k}.json")) for k in ("kpgp3_self", "kpgp6_on3", "kp1_on3", "kpmtx_on3", "kpgp3_ont", "kpgp6_ont", "kp1_ont", "kpmtx_ont")}
def series(k):
    b = np.array(prof[k]["bins"]); m = (b[:, 0] >= X0 - 100) & (b[:, 0] <= X1)
    return b[m, 0] + 50, b[m, 2]
# narrow mouse-unmappable strips, detected from the Illumina bins (display only): both mouse < 0.15, both GP >= 0.30,
# outside MOUSE_ZERO; bins joined across at most one 100-bp gap; strips >= 100 bp kept
STRIP_RULE = "kp1_on3 < 0.15 & kpmtx_on3 < 0.15 & kpgp3_self >= 0.30 & kpgp6_on3 >= 0.30; outside MOUSE_ZERO; <=1-bin gap merged; width >= 100 bp"
def mouse_strips():
    B = {k: {int(s): v for s, _, v in prof[k]["bins"]} for k in ("kp1_on3", "kpmtx_on3", "kpgp3_self", "kpgp6_on3")}
    hit = [s for s in sorted(B["kp1_on3"]) if X0 <= s and s + 99 <= X1 and not (MOUSE_ZERO[0] <= s <= MOUSE_ZERO[1])
           and B["kp1_on3"][s] < 0.15 and B["kpmtx_on3"][s] < 0.15 and B["kpgp3_self"][s] >= 0.30 and B["kpgp6_on3"][s] >= 0.30]
    out = []
    for s in hit:
        if out and s - out[-1][1] <= 101: out[-1][1] = s + 99   # adjacent (1) or one empty bin between (101)
        else: out.append([s, s + 99])
    return [tuple(x) for x in out if x[1] - x[0] + 1 >= 100]
# B (Illumina): every detected strip. B2 (ONT): only strips where mouse ONT depth also drops (both mouse ONT < 0.15 in >= 1 bin),
# i.e. mouse lacks the sequence. A strip without the ONT drop (93.5 kb, inside the repeat copy) is short-read mapping ambiguity
# (MAPQ < 20), not Kp-GP-specific sequence; the panel B note says so. Checked in verify_strips.py.
def ont_confirmed(a, b):
    O = {k: {int(s): v for s, _, v in prof[k]["bins"]} for k in ("kp1_ont", "kpmtx_ont")}
    return any(O["kp1_ont"][s] < 0.15 and O["kpmtx_ont"][s] < 0.15 for s in range(a, b + 1, 100))
STRIPS_ALL = mouse_strips()
STRIPS = [x for x in STRIPS_ALL if ont_confirmed(*x)]          # drawn in B and B2
STRIPS_B_ONLY = [x for x in STRIPS_ALL if x not in STRIPS]      # drawn in B only
print("strips (B and B2):", ", ".join(f"{a:,}-{b:,}" for a, b in STRIPS), "| B only:", ", ".join(f"{a:,}-{b:,}" for a, b in STRIPS_B_ONLY))
genes =[r for r in csv.DictReader(open(f"{CNV}/ahp1_region_genes.tsv"), delimiter="\t") if r["gene_id"] in GENE_NAMES]
ev = {k: json.load(open(f"{CNV}/evidence_{k}.json")) for k in ("kpgp3_self", "kpgp6_on3", "kp1_on3", "kpmtx_on3")}
def pairs(k):
    return sorted({(p, q) for p, c, q in ev[k]["disc_targets"] if c == "scaffold_18" and 78500 < p < 80000 and 93000 < q < 95000})
ont = {s: json.load(open(f"{HERE}/ont_reads_{s}.json"))["reads"] for s in ("kpgp3", "kpgp6")}

MM = 1 / 25.4
FIG_W, FIG_H = 174.0, 186.0                              # declared final print size in mm
fig = plt.figure(figsize=(FIG_W * MM, FIG_H * MM))
gs = fig.add_gridspec(5, 1, height_ratios=[1.25, 3.45, 1.15, 2.7, 1.5], hspace=0.30,
                      left=22.0 / FIG_W, right=1 - 1.5 / FIG_W, top=1 - 5.0 / FIG_H, bottom=10.0 / FIG_H)
axA, axB, axB2, axC, axD = [fig.add_subplot(gs[i]) for i in range(5)]
for ax in (axA, axB, axB2, axC, axD):
    ax.set_xlim(X0, X1)
    for d in DEL: ax.axvline(d, color=ACC, lw=0.7, ls=(0, (3, 2)), zorder=1)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
for ax in (axA, axB, axB2, axC): ax.tick_params(labelbottom=False)
for ax in (axA, axC, axD): ax.spines["left"].set_visible(False); ax.set_yticks([])
for ax in (axA, axC): ax.spines["bottom"].set_visible(False); ax.tick_params(bottom=False)

# ---------------- A: genes ----------------
axA.set_ylim(-1.6, 2.2)
axA.plot([X0, X1], [0, 0], color="#999", lw=0.8, zorder=0)
for i, g in enumerate(sorted(genes, key=lambda r: int(r["start"]))):
    s, e, st = int(g["start"]), int(g["end"]), g["strand"]
    inside = g["position_vs_deletion"] == "INSIDE"
    fc, ec = (ACC_F, ACC) if inside else ("#d9d9d9", "#999")
    if g["type"] == "tRNA":
        axA.add_patch(Rectangle((s, -0.28), max(e - s, 120), 0.56, fc=fc, ec=ec, lw=0.8))
    else:
        h, head = 0.32, min(400, (e - s) * 0.35)
        if st == "+": pts = [(s, -h), (e - head, -h), (e - head, -h * 1.6), (e, 0), (e - head, h * 1.6), (e - head, h), (s, h)]
        else: pts = [(e, -h), (s + head, -h), (s + head, -h * 1.6), (s, 0), (s + head, h * 1.6), (s + head, h), (e, h)]
        pts = [(max(X0, x), y) for x, y in pts]
        axA.add_patch(Polygon(pts, closed=True, fc=fc, ec=ec, lw=0.8, zorder=2))
    lab = GENE_NAMES[g["gene_id"]]
    y = 0.75 if i % 2 == 0 else -0.95
    axA.text((max(s, X0) + e) / 2, y, lab, ha="center", va="center", fontsize=7.5,
             fontweight="bold" if g["gene_id"] in BOLD else "normal", color=ACC if inside else "#666", fontstyle="normal")
for a, b in REP:
    axA.add_patch(Polygon([(a, -1.45), (b, -1.45), ((a + b) / 2, -1.15)], closed=True, fc="#222", ec="#222", lw=0.5))
axA.text(REP[0][0] - 150, -1.3, "316-bp repeat", ha="right", va="center", fontsize=7.5, color="#666")
axA.plot([DEL[0], DEL[0], DEL[1], DEL[1]], [1.55, 1.75, 1.75, 1.55], color=ACC, lw=1.0)
axA.text((DEL[0] + DEL[1]) / 2, 1.95, "13.8-kb deletion", ha="center", va="bottom", fontsize=8.0, color=ACC, fontweight="bold")
axA.text(-0.105, 1.0, "A", transform=axA.transAxes, fontsize=9, fontweight="bold", va="top", ha="left")

# ---------------- B: Illumina depth ----------------
axB.add_patch(Rectangle((MOUSE_ZERO[0], 0), MOUSE_ZERO[1] - MOUSE_ZERO[0], 2.05, fc=BAND, ec="none", zorder=0))
for a, b in STRIPS_ALL: axB.add_patch(Rectangle((a, 0), b - a, 2.05, fc=BAND, ec="none", zorder=0))
axB.text(X1, 1.98, "grey = mouse reads not mapped", ha="right", va="top",
         fontsize=7.5, color="#4b5563")
for k, s, col, ls in (("kp1_on3", "kp1", MOUSE_C, "-"), ("kpmtx_on3", "kpmtx", MOUSE_C, (0, (2.5, 1.5))),
                      ("kpgp3_self", "kpgp3", GP_C, "-"), ("kpgp6_on3", "kpgp6", GP_C, (0, (2.5, 1.5)))):
    x, y = series(k); axB.plot(x, y, color=col, ls=ls, lw=0.75, label=STRAIN_NAMES[s])
axB.axhline(1.0, color="#999", lw=0.6)
axB.set_ylim(0, 2.05); axB.set_yticks([0, 0.5, 1.0, 1.5, 2.0]); axB.set_ylabel("Illumina depth\n(self median = 1)")
axB.legend(fontsize=7.5, frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0), handlelength=2.2, columnspacing=1.6)
inner = lambda k: float(np.mean([v for x, v in zip(*series(k)) if 84000 <= x - 50 < 93000]))
txt = "84–93 kb mean:  " + "   ".join(f"{STRAIN_NAMES[s]} {inner(k):.2f}" for k, s in (("kpgp3_self", "kpgp3"), ("kpgp6_on3", "kpgp6"), ("kp1_on3", "kp1"), ("kpmtx_on3", "kpmtx")))
axB.text(-0.105, 1.09, "B", transform=axB.transAxes, fontsize=9, fontweight="bold", va="top", ha="left")

# ---------------- B2: ONT depth ----------------
axB2.add_patch(Rectangle((MOUSE_ZERO[0], 0), MOUSE_ZERO[1] - MOUSE_ZERO[0], 2.0, fc=BAND, ec="none", zorder=0))
for a, b in STRIPS: axB2.add_patch(Rectangle((a, 0), b - a, 2.0, fc=BAND, ec="none", zorder=0))
for k, s, col, ls in (("kp1_ont", "kp1", MOUSE_C, "-"), ("kpmtx_ont", "kpmtx", MOUSE_C, (0, (2.5, 1.5))),
                      ("kpgp3_ont", "kpgp3", GP_C, "-"), ("kpgp6_ont", "kpgp6", GP_C, (0, (2.5, 1.5)))):
    x, y = series(k); axB2.plot(x, y, color=col, ls=ls, lw=0.6)
axB2.axhline(1.0, color="#999", lw=0.6); axB2.set_ylim(0, 2.0); axB2.set_yticks([0, 1, 2])
axB2.set_ylabel("ONT depth\n(self median = 1)", fontsize=7.5)

# ---------------- C: ONT deletion reads ----------------
rows = [(s, r) for s in ("kpgp3", "kpgp6") for r in sorted(ont[s], key=lambda r: (r["gap"][1] - r["gap"][0] < 10000, r["segments"][0][0]))]
for i, (s, r) in enumerate(rows):
    y = -i
    junction = r["gap"][1] - r["gap"][0] > 10000
    col = GP_C if junction else "#999"
    (a0, a1), (b0, b1) = r["segments"]
    axC.plot([max(a0, X0), a1], [y, y], color=col, lw=1.1 if s == "kpgp3" else 1.1, solid_capstyle="butt", alpha=1 if s == "kpgp3" else 0.6)
    axC.plot([b0, min(b1, X1)], [y, y], color=col, lw=1.1, solid_capstyle="butt", alpha=1 if s == "kpgp3" else 0.6)
    axC.plot([a1, b0], [y, y], color=col, lw=0.4, ls=(0, (1, 1.5)), alpha=0.8)
n3 = len(ont["kpgp3"]); n6 = len(ont["kpgp6"])
j3 = sum(1 for r in ont["kpgp3"] if r["gap"][1] - r["gap"][0] > 10000)
axC.set_ylim(-len(rows) - 1.6, 6.1)
axC.text(X0, 5.9, "ONT reads crossing the deletion · blue = junction, grey = other\nsolid = aligned, dotted = reference skipped",
         ha="left", va="top", fontsize=7.5, linespacing=1.3)
axC.text(X1, 5.9, f"{STRAIN_NAMES['kpgp3']}  {j3} junction + {n3 - j3} other", ha="right", va="top", fontsize=7.5, color=GP_C)
axC.text(X1, -len(rows) - 0.1, f"{STRAIN_NAMES['kpgp6']}  {n6} junction", ha="right", va="top", fontsize=7.5, color=GP_C, alpha=0.55)
axC.text(-0.105, 0.99, "C", transform=axC.transAxes, fontsize=9, fontweight="bold", va="top", ha="left")

# ---------------- D: Illumina discordant pairs ----------------
axD.set_ylim(-4.3, 1.6)
for j, (k, s) in enumerate((("kpgp3_self", "kpgp3"), ("kpgp6_on3", "kpgp6"), ("kp1_on3", "kp1"), ("kpmtx_on3", "kpmtx"))):
    P = pairs(k); base = [0.0, -1.2, -2.4, -3.6][j]
    col = GP_C if s.startswith("kpgp") else MOUSE_C
    for p, q in P:
        axD.add_patch(Arc(((p + q) / 2, base), q - p, 2.2 if j < 2 else 1, theta1=0, theta2=180, color=col, lw=0.5, alpha=0.8))
    axD.text(X0 + 100, base + 0.05, f"{STRAIN_NAMES[s]}  {len(P)} pairs", fontsize=7.5, color=col, va="bottom", ha="left")
axD.text(X0, 1.2, "Illumina FR pairs spanning the junction, one arc per pair", ha="left", va="bottom", fontsize=7.5, color="#222")
axD.set_xticks(range(76000, 97001, 2000)); axD.set_xticklabels([f"{t / 1000:.0f}" for t in range(76000, 97001, 2000)])
axD.set_xlabel(f"{STRAIN_NAMES['kpgp3']} scaffold_18 position (kb)")
axD.text(-0.105, 1.0, "D", transform=axD.transAxes, fontsize=9, fontweight="bold", va="top", ha="left")
fig.savefig(f"{HERE}/hog1_ahp1_deletion_figure.svg")      # no bbox_inches: the canvas above is the printed size
fig.savefig(f"{HERE}/hog1_ahp1_deletion_figure.png", dpi=600)

# ---------------- data for reproduction ----------------
with open(f"{HERE}/deletion_figure_data.tsv", "w") as fh:
    fh.write("# section 1: 100-bp bin normalized depth (bin centre), from ahp1_cnv/profile_*.json (values copied, not recomputed)\n")
    keys = ["kpgp3_self", "kpgp6_on3", "kp1_on3", "kpmtx_on3", "kpgp3_ont", "kpgp6_ont", "kp1_ont", "kpmtx_ont"]
    fh.write("bin_start\t" + "\t".join(f"{k}(median {prof[k]['median']:g})" for k in keys) + "\n")
    b0 = np.array(prof[keys[0]]["bins"])
    for i, st in enumerate(b0[:, 0]):
        if X0 - 100 <= st <= X1:
            fh.write(f"{int(st)}\t" + "\t".join(f"{prof[k]['bins'][i][2]:.4f}" for k in keys) + "\n")
    fh.write("\n# section 2: genes (ahp1_cnv/ahp1_region_genes.tsv)\ngene_id\tlabel\ttype\tstart\tend\tstrand\tposition_vs_deletion\n")
    for g in genes: fh.write(f"{g['gene_id']}\t{GENE_NAMES[g['gene_id']]}\t{g['type']}\t{g['start']}\t{g['end']}\t{g['strand']}\t{g['position_vs_deletion']}\n")
    fh.write(f"\n# section 3: features\ndeleted_interval\t{DEL[0]}\t{DEL[1]}\nrepeat_copy_1\t{REP[0][0]}\t{REP[0][1]}\nrepeat_copy_2\t{REP[1][0]}\t{REP[1][1]}\n"
             f"not_in_mouse_assemblies_blastn\t{MOUSE_ABSENT[0]}\t{MOUSE_ABSENT[1]}\nmouse_depth_near_zero_band\t{MOUSE_ZERO[0]}\t{MOUSE_ZERO[1]}\n"
             + f"# mouse-unmappable narrow strips (grey, panels B and B2), detected from the Illumina bins: {STRIP_RULE}\n"
             + "# drawn in B and B2 if both mouse ONT profiles are also < 0.15 in >= 1 bin of the strip; otherwise B only (see verify_strips.py)\n"
             + "".join(f"mouse_unmappable_strip\t{a}\t{b}\t{'drawn_B_and_B2' if (a, b) in STRIPS else 'drawn_B_only_mouse_ONT_depth_normal'}\n" for a, b in STRIPS_ALL))
    fh.write("\n# section 4: ONT deletion reads (ont_reads_*.json)\nstrain\tread\ttype\tseg1_start\tseg1_end\tseg2_start\tseg2_end\tgap_bp\n")
    for s in ("kpgp3", "kpgp6"):
        for r in ont[s]:
            (a0, a1), (b0_, b1) = r["segments"]; fh.write(f"{s}\t{r['read']}\t{r['type']}\t{a0}\t{a1}\t{b0_}\t{b1}\t{r['gap'][1] - r['gap'][0]}\n")
    fh.write("\n# section 5: Illumina discordant FR pairs (unique read/mate positions; evidence_*.json)\nstrain\tread_pos\tmate_pos\tspan\n")
    for k, s in (("kpgp3_self", "kpgp3"), ("kpgp6_on3", "kpgp6"), ("kp1_on3", "kp1"), ("kpmtx_on3", "kpmtx")):
        P = pairs(k)
        if not P: fh.write(f"{s}\tnone\tnone\t0\n")
        for p, q in P: fh.write(f"{s}\t{p}\t{q}\t{q - p}\n")
    fh.write("\n# section 6: interior (84-93 kb) mean of the plotted bins\n" + "\n".join(f"{k}\t{inner(k):.4f}" for k in keys) + "\n")
print("saved")
