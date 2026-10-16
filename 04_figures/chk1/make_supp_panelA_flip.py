#!/usr/bin/env python3
"""chk1_deletion.{svg,png} - the 1.9-kb CHK1 deletion in the guinea-pig clade, on the Kp-M (kp1) assembly.
Built in the same visual grammar as ahp1_cnv/figure/make_figure.py (HOG1-AHP1): the same five panels in the same
order (A genes [+ an unlettered protein row], B Illumina depth, an unlettered ONT depth row, C individual ONT
reads, D Illumina junction pairs), the same strain order, colours, line styles and depth axis, the deletion
boundaries as dashed accent verticals on every panel, a declared print size and no bbox_inches crop.
Only existing results are read (no re-estimation):
  depth      chk1/work/q_{kp1,kpmtx,kpgp3,kpgp6}_{il,ont}.json (100-bp bins, MAPQ>=20, / own genome median)
  genes      chk1/work/genes_labels.tsv (kp1 funannotate; labels = reciprocal best blastp hit in S. cerevisiae R64)
  ONT reads  the ont_reads field of the same q_*_ont.json (classes junction/intact/other, chk1/reconcile_ont.py)
  pairs      the pair_list field of the same q_*_il.json (the dedicated chk1 definition: fixed flanking windows
             417.5-420.99 kb and 422.86-426 kb) - 83 / 149 / 0 / 0, NOT the sv_scan definition (90 / 164)
The interval drawn is the validated deletion boundary 420,988-422,860 (1,873 bp), never the scan window.

In the protein row of panel A the residue axis runs 570 at the left to 1 at the right: CHK1 is on the minus
strand, so residue 1 lies at the high genomic coordinate and residue 570 at the low one, and the row then
points the same way as the gene map above it. N and C are written at the two ends, N at the right and C at
the left. The ONT read sample is seeded (random.Random(1))."""
import json, csv, os, random
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, Arc
plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "font.family": "Arial", "font.size": 7.5,
                     "axes.edgecolor": "#666", "xtick.color": "#666", "ytick.color": "#666",
                     "axes.labelcolor": "#222", "text.color": "#222"})
HERE = os.path.dirname(os.path.abspath(__file__)); W = f"{HERE}/../work"

# ---------------- naming and style (identical to the HOG1-AHP1 figure) ----------------
STRAIN_NAMES = {"kp1": "Kp-M", "kpmtx": "Kp-MTX", "kpgp3": "Kp-GP", "kpgp6": "Kp-GP2"}
MOUSE_C, GP_C, ACC, ACC_F = "#F4A261", "#3a7ebf", "#A32E2E", "#FDECEC"
BAND = "#e5e7eb"
ST = [("kp1", MOUSE_C, "-"), ("kpmtx", MOUSE_C, (0, (2.5, 1.5))), ("kpgp3", GP_C, "-"), ("kpgp6", GP_C, (0, (2.5, 1.5)))]
X0, X1 = 414000, 430000
DEL = (420988, 422860)                   # validated deletion boundary, 1,873 bp (kp1 coordinates)
PROT_LEN, PKIN, REMOVED = 570, (24, 297), (1, 520)
BOLD = {"CHK1"}
N_ONT = 5                                # ONT reads drawn per strain

q = lambda s, t: json.load(open(f"{W}/q_{s}_{t}.json"))
def series(s, t):
    r = q(s, t); x = (r["profile_start"] + 50 + 100 * np.arange(len(r["profile100"]))).astype(float)
    y = np.array(r["profile100"], dtype=float); m = (x >= X0) & (x <= X1)
    return x[m], y[m]

# ---------------- canvas: declared print size, in mm ----------------
MM = 1 / 25.4
FIG_W, FIG_H = 174.0, 195.0
LEFT, RIGHT = 22.0, 4.0                                  # the 430 tick label is centred on the right edge
fig = plt.figure(figsize=(FIG_W * MM, FIG_H * MM))
def box(top_mm, h_mm, inset_mm=0.0):                     # top-down mm -> figure fraction
    return fig.add_axes([(LEFT + inset_mm) / FIG_W, 1 - (top_mm + h_mm) / FIG_H,
                         (FIG_W - LEFT - RIGHT - 2 * inset_mm) / FIG_W, h_mm / FIG_H])
axA  = box(5.0, 20.0)
axP  = box(29.0, 8.0, inset_mm=20.0)                     # inset and framed: NOT on the genomic axis above
axB  = box(48.0, 42.0)
axB2 = box(96.0, 17.0)
axC  = box(117.0, 46.0)
axD  = box(169.0, 14.0)
for ax in (axA, axB, axB2, axC, axD):
    ax.set_xlim(X0, X1)
    for d in DEL: ax.axvline(d, color=ACC, lw=0.7, ls=(0, (3, 2)), zorder=1)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
for ax in (axA, axB, axB2, axC): ax.tick_params(labelbottom=False)
for ax in (axA, axC, axD): ax.spines["left"].set_visible(False); ax.set_yticks([])
for ax in (axA, axC): ax.spines["bottom"].set_visible(False); ax.tick_params(bottom=False)
L = lambda ax, s, y=1.0: ax.text(-0.105, y, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top", ha="left")

# ---------------- A: genes ----------------
genes = [g for g in csv.DictReader(open(f"{W}/genes_labels.tsv"), delimiter="\t")
         if int(g["end"]) >= X0 and int(g["start"]) <= X1]
axA.set_ylim(-1.6, 2.2)
axA.plot([X0, X1], [0, 0], color="#999", lw=0.8, zorder=0)
for i, g in enumerate(sorted(genes, key=lambda r: int(r["start"]))):
    s, e, st = int(g["start"]), int(g["end"]), g["strand"]
    inside = not (e < DEL[0] or s > DEL[1])
    fc, ec = (ACC_F, ACC) if inside else ("#d9d9d9", "#999")
    h, head = 0.32, min(300, (e - s) * 0.35)
    if st == "+": pts = [(s, -h), (e - head, -h), (e - head, -h * 1.6), (e, 0), (e - head, h * 1.6), (e - head, h), (s, h)]
    else: pts = [(e, -h), (s + head, -h), (s + head, -h * 1.6), (s, 0), (s + head, h * 1.6), (s + head, h), (e, h)]
    pts = [(min(X1, max(X0, x)), y) for x, y in pts]
    axA.add_patch(Polygon(pts, closed=True, fc=fc, ec=ec, lw=0.8, zorder=2))
    lab = g["label"]
    if lab:
        axA.text((min(X1, max(s, X0)) + min(e, X1)) / 2, 0.75 if i % 2 == 0 else -0.95, lab, ha="center", va="center",
                 fontsize=7.5, fontweight="bold" if lab in BOLD else "normal", color=ACC if inside else "#666")
axA.plot([DEL[0], DEL[0], DEL[1], DEL[1]], [1.55, 1.75, 1.75, 1.55], color=ACC, lw=1.0)
axA.text((DEL[0] + DEL[1]) / 2, 1.95, "1.9-kb deletion", ha="center", va="bottom", fontsize=8.0, color=ACC, fontweight="bold")
L(axA, "A")

# ---------------- protein row (continuation of A, no letter) ----------------
# residue axis reversed: CHK1 is on the minus strand, so residue 1 sits at the high genomic coordinate (right)
axP.set_xlim(PROT_LEN + 14, -14); axP.set_ylim(-0.15, 1.65)
axP.set_facecolor("#fafafa")
for sp in axP.spines.values(): sp.set_visible(True); sp.set_edgecolor("#aaa"); sp.set_linewidth(0.5)
axP.add_patch(Rectangle((1, 0.52), PROT_LEN - 1, 0.32, fc="#f3f4f6", ec="#999", lw=0.6))
axP.add_patch(Rectangle((REMOVED[0], 0.52), REMOVED[1] - REMOVED[0], 0.32, fc=ACC_F, ec=ACC, lw=0.8))
axP.add_patch(Rectangle((PKIN[0], 0.58), PKIN[1] - PKIN[0], 0.20, fc="#9ca3af", ec="none"))
axP.text((PKIN[0] + PKIN[1]) / 2, 0.94, "Pkinase (PF00069)", ha="center", va="bottom", fontsize=7.5, color="#4b5563")
axP.text((REMOVED[0] + REMOVED[1]) / 2, 0.42, f"residues {REMOVED[0]}–{REMOVED[1]} removed", ha="center", va="top", fontsize=7.5, color=ACC)
axP.set_xticks([1, 100, 200, 300, 400, 500, PROT_LEN]); axP.tick_params(axis="x", length=2, labelsize=7.5, pad=2)
axP.text(-7.5, 0.68, "N", ha="center", va="center", fontsize=7.5, color="#4b5563", fontweight="bold")
axP.text(PROT_LEN + 7.5, 0.68, "C", ha="center", va="center", fontsize=7.5, color="#4b5563", fontweight="bold")
axP.set_yticks([])
axP.text(-0.03, 0.5, f"CHK1 protein ({PROT_LEN} aa)", transform=axP.transAxes,
         ha="right", va="center", fontsize=7.5, color="#4b5563", linespacing=1.25)

# ---------------- B: Illumina depth ----------------
lo = min(q(s, "il")["noise"]["q"]["5"] for s, _, _ in ST); hi = max(q(s, "il")["noise"]["q"]["95"] for s, _, _ in ST)
axB.axhspan(lo, hi, color=BAND, lw=0, zorder=0)
for s, col, ls in ST:
    x, y = series(s, "il"); axB.plot(x, y, color=col, ls=ls, lw=0.75, label=STRAIN_NAMES[s])
axB.axhline(1.0, color="#999", lw=0.6)
axB.set_ylim(-0.09, 2.05); axB.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
axB.set_ylabel("Illumina depth\n(self median = 1)")
axB.legend(fontsize=7.5, frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0), handlelength=2.2, columnspacing=1.6)
axB.text(X1, 1.98, "grey = genome-wide 5–95%", ha="right", va="top", fontsize=7.5, color="#4b5563")
for _ax in (axB, axB2):                                  # a guinea-pig depth of exactly 0 is a measured zero, not a gap
    _ax.plot([(DEL[0] + DEL[1]) / 2], [0], marker="_", ms=13, mew=1.5, color=GP_C, zorder=4)
L(axB, "B", 1.09)

# ---------------- ONT depth (continuation of B, no letter) ----------------
for s, col, ls in ST:
    x, y = series(s, "ont"); axB2.plot(x, y, color=col, ls=ls, lw=0.6)
axB2.axhline(1.0, color="#999", lw=0.6); axB2.set_ylim(-0.09, 2.0); axB2.set_yticks([0, 1, 2])
axB2.set_ylabel("ONT depth\n(self median = 1)")

# ---------------- C: individual ONT reads ----------------
rng = random.Random(1); yy = 0.0; CNT = {}
axC.text(X0, 3.1, f"ONT reads at the locus (≤ {N_ONT} per strain) · colour = strain · solid = aligned, dotted = reference skipped",
         ha="left", va="bottom", fontsize=7.5)
for s, col, ls in ST:
    r = q(s, "ont"); c = r["ont_counts"]; CNT[s] = c
    assert set(x["cls"] for x in r["ont_reads"]) <= {"junction", "intact", "other"}, "ONT class names changed"
    dele = [x for x in r["ont_reads"] if x["cls"] == "junction"]
    intact = [x for x in r["ont_reads"] if x["cls"] == "intact"]
    pick = rng.sample(dele, min(N_ONT, len(dele))) + rng.sample(intact, min(max(0, N_ONT - len(dele)), len(intact)))
    assert len(dele) == c["junction"] and len(intact) == c["intact"], "counts and drawn reads disagree"
    axC.text(X0, yy + 0.6, f"{STRAIN_NAMES[s]}   {c['junction']} junction, {c['intact']} intact "
             f"/ {c['junction'] + c['intact'] + c['other']}", ha="left", va="bottom", fontsize=7.5, color=col)
    yy -= 1.0
    if not pick:
        axC.text((DEL[0] + DEL[1]) / 2, yy, "no ONT read crosses the interval", ha="center", va="center", fontsize=7.5, color="#666"); yy -= 1
    for x in sorted(pick, key=lambda z: min(b[0] for b in z["blocks"])):
        bl = sorted(x["blocks"])
        for a, b in bl: axC.plot([max(a, X0), min(b, X1)], [yy, yy], color=col, lw=1.1, solid_capstyle="butt")
        for (a1, b1), (a2, b2) in zip(bl, bl[1:]):
            if a2 > b1: axC.plot([b1, a2], [yy, yy], color=col, lw=0.4, ls=(0, (1, 1.5)))
        yy -= 1
    yy -= 3.0
axC.set_ylim(yy + 1.6, 4.1)
L(axC, "C")

# ---------------- D: Illumina junction pairs ----------------
axD.set_ylim(-4.6, 1.3)
for j, (s, col, ls) in enumerate(ST):
    P = sorted({tuple(p) for p in q(s, "il")["pair_list"]}); base = [0.0, -1.3, -2.6, -3.9][j]
    col = GP_C if s.startswith("kpgp") else MOUSE_C
    for p, r_ in P:
        axD.add_patch(Arc(((p + r_) / 2, base), r_ - p, 0.9, theta1=0, theta2=180, color=col, lw=0.4, alpha=0.5))
    axD.text(X0 + 150, base + 0.05, f"{STRAIN_NAMES[s]}  {len(P)} pairs", fontsize=7.5, color=col, va="bottom", ha="left")
axD.text(X0, 1.4, "Illumina FR pairs spanning the junction, one arc per pair", ha="left", va="bottom", fontsize=7.5, color="#222")
axD.set_xticks(range(414000, 430001, 2000)); axD.set_xticklabels([f"{t / 1000:.0f}" for t in range(414000, 430001, 2000)])
axD.set_xlabel(f"{STRAIN_NAMES['kp1']} (kp1) scaffold_6 position (kb)")
L(axD, "D")

OUT = "chk1_deletion_panelA_flip"
fig.savefig(f"{HERE}/{OUT}.svg")                   # no bbox_inches: the canvas above is the printed size
fig.savefig(f"{HERE}/{OUT}.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})

# ---------------- data for reproduction ----------------
with open(os.devnull, "w") as fh:          # the existing chk1_deletion_data.tsv is NOT rewritten by this script
    fh.write("# section 1: 100-bp bin normalized depth (bin centre), from chk1/work/q_*.json (values copied, not recomputed)\n")
    keys = [(s, t) for t in ("il", "ont") for s, _, _ in ST]
    fh.write("bin_start\t" + "\t".join(f"{s}_{t}(median {q(s, t)['median']:g})" for s, t in keys) + "\n")
    xs, _ = series(*keys[0])
    cols = {k: dict(zip(*series(*k))) for k in keys}
    for x in xs: fh.write(f"{int(x - 50)}\t" + "\t".join(f"{cols[k].get(x, float('nan')):.4f}" for k in keys) + "\n")
    fh.write("\n# section 2: genes (chk1/work/genes_labels.tsv)\nstart\tend\tstrand\tkp1_gene\tlabel\tinside_deletion\n")
    for g in sorted(genes, key=lambda r: int(r["start"])):
        ins = not (int(g["end"]) < DEL[0] or int(g["start"]) > DEL[1])
        fh.write(f"{g['start']}\t{g['end']}\t{g['strand']}\t{g['kp1_gene']}\t{g['label'] or '-'}\t{ins}\n")
    fh.write(f"\n# section 3: features (kp1 coordinates)\nvalidated_deletion_boundary\t{DEL[0]}\t{DEL[1]}\t{DEL[1] - DEL[0] + 1} bp\n"
             f"chk1_protein_length_aa\t{PROT_LEN}\npkinase_PF00069_aa\t{PKIN[0]}\t{PKIN[1]}\nresidues_removed_aa\t{REMOVED[0]}\t{REMOVED[1]}\n")
    fh.write("\n# section 4: ONT read counts (q_*_ont.json field ont_counts; 'other' = does not reach both sides)\nstrain\tjunction\tintact\tother\n")
    for s, _, _ in ST: fh.write(f"{s}\t{CNT[s]['junction']}\t{CNT[s]['intact']}\t{CNT[s]['other']}\n")
    fh.write("\n# section 5: Illumina junction pairs (dedicated chk1 definition; q_*_il.json field pair_list)\nstrain\tread_pos\tmate_pos\tspan\n")
    for s, _, _ in ST:
        P = sorted({tuple(p) for p in q(s, "il")["pair_list"]})
        if not P: fh.write(f"{s}\tnone\tnone\t0\n")
        for p, r_ in P: fh.write(f"{s}\t{p}\t{r_}\t{r_ - p}\n")
print("saved")
