#!/usr/bin/env python3
"""Promoter cis-element changes (STRE / YAP1) at stress genes, mouse clade vs guinea-pig clade.
Reads stre_turnover_events.tsv only. Declared print size, no bbox_inches crop, so a matplotlib font size of n is
n pt on the page. Two visual axes, both in the legend: shape says whether a site differs between the clades
(triangle / diamond) or is shared by all nine strains (tick); fill says whether the site is present or absent."""
import csv, os, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.family": "Arial", "svg.fonttype": "none", "font.size": 7.5, "text.color": "#222",
                     "axes.edgecolor": "#666", "xtick.color": "#666", "axes.labelcolor": "#222"})
MOUSE_C, GP_C = "#F4A261", "#3a7ebf"
# No gene list: every locus of the 70-gene stress set that meets the criterion is drawn, and the rows are in
# alphabetical order. AHP1 carries an asterisk because its promoter has two guinea-pig copies rather than three,
# which the caption explains; it is not boxed or shaded, because a highlight would read as importance and the
# point of the figure is that these eleven are a result of the criterion, not a selection.
MARK = {"STRE": "^", "YAP1": "D"}
MS = {"STRE": 5.0, "YAP1": 4.2}
MM = 1 / 25.4
FIG_W, FIG_H = 174.0, 171.0                               # declared final print size in mm

rows = list(csv.DictReader([l for l in open(f"{HERE}/stre_turnover_events.tsv") if not l.startswith("#")], delimiter="\t"))
centre = lambda p: sum(int(x) for x in p.split("..")) / 2
N_SHARED = sum(r["event"] == "shared" for r in rows)
N_DIFF = sum(r["event"] != "shared" for r in rows)
assert N_SHARED + N_DIFF == len(rows), (N_SHARED, N_DIFF, len(rows))
assert N_DIFF == 17, f"expected 17 clade differences, got {N_DIFF}"
assert N_SHARED == 23, f"expected 23 shared sites, got {N_SHARED}"
GENES = sorted({r["gene"] for r in rows})
assert len(GENES) == 11, f"expected 11 genes, got {len(GENES)}"
_PER = collections.Counter(r["gene"] for r in rows if r["event"] != "shared")
DIST = collections.Counter(_PER.values())                  # events per gene -> how many genes
assert sum(k * v for k, v in DIST.items()) == N_DIFF, (DIST, N_DIFF)
assert sum(DIST.values()) == len(GENES), (DIST, len(GENES))
ORDER = GENES                                              # alphabetical, no exception: the rows are a result,
                                                           # not a ranking, so no row is given a position of its own
drawn = {"tick": 0, "symbol": 0}      # counted while drawing and checked against the table at the end

fig = plt.figure(figsize=(FIG_W * MM, FIG_H * MM))
ax = fig.add_axes([14.0 / FIG_W, 11.0 / FIG_H, (FIG_W - 16.0) / FIG_W, (FIG_H - 11.0 - 13.0) / FIG_H])
ROW, DY, LAB = 1.75, 0.34, 0.34                           # row pitch; half-gap between the two tracks; label drop
y = 0.0
if True:
    for k, g in enumerate(ORDER):
        y -= 1.0 if k == 0 else ROW
        ym, yg = y + DY, y - DY
        ax.text(-712, y, g + ("*" if g == "AHP1" else ""), ha="right", va="center", fontsize=8.0, fontweight="bold")
        ax.plot([-700, 0], [ym, ym], color=MOUSE_C, lw=1.0, solid_capstyle="butt", zorder=1)
        ax.plot([-700, 0], [yg, yg], color=GP_C, lw=1.0, solid_capstyle="butt", zorder=1)
        diffs = sorted((r for r in rows if r["gene"] == g and r["event"] != "shared"), key=lambda r: centre(r["pos_from_ATG"]))
        drop = {}                                          # labels closer than 140 bp get a second row
        prev_x, prev_lvl = -1e9, 1
        for r in diffs:
            cx = centre(r["pos_from_ATG"])
            lvl = 2 if cx - prev_x < 140 and prev_lvl == 1 else 1
            drop[id(r)] = lvl; prev_x, prev_lvl = cx, lvl
        for r in (r for r in rows if r["gene"] == g):
            x, m, ev = centre(r["pos_from_ATG"]), r["motif"], r["event"]
            if ev == "shared":                             # a tick, not a symbol: no clade difference to read
                for yy, c in ((ym, MOUSE_C), (yg, GP_C)):
                    ax.plot([x, x], [yy - 0.13, yy + 0.13], color=c, lw=0.9, solid_capstyle="butt", zorder=2)
                    drawn["tick"] += 1
                continue
            mouse_on, gp_on = ev == "loss", ev == "gain"
            ax.plot([x, x], [yg, ym], color="#555", lw=0.5, ls=(0, (1.5, 1.2)), zorder=1.5)
            for yy, c, on in ((ym, MOUSE_C, mouse_on), (yg, GP_C, gp_on)):
                ax.plot(x, yy, MARK[m], ms=MS[m] + 1.0, mfc=c if on else "white", mec=c, mew=1.1, zorder=3)
                drawn["symbol"] += 1
            lab = f"{m} {ev} {r['pos_from_ATG'].split('..')[0].replace('-', '−')}"
            if "indel" in r["cause"]:
                lab += " (indel)"
            ly = yg - LAB - 0.42 * (drop[id(r)] - 1)
            ax.text(x, ly, lab, ha="center", va="top", fontsize=7.5, color="#222",
                    fontweight="bold" if ev == "gain" else "normal")
assert drawn["tick"] == 2 * N_SHARED, (drawn, N_SHARED)     # one tick per track
assert drawn["symbol"] == 2 * N_DIFF, (drawn, N_DIFF)       # one symbol per track
y_end = y - 1.9

ax.set_xlim(-860, 25); ax.set_ylim(y_end, 1.5)
ax.set_yticks([])
for sp in ("left", "right", "top"):
    ax.spines[sp].set_visible(False)
ax.spines["bottom"].set_bounds(-700, 0)
ax.set_xticks(range(-700, 1, 100))
ax.set_xticklabels([("ATG" if t == 0 else f"−{-t}") for t in range(-700, 1, 100)], fontsize=7.5)
ax.set_xlabel("position relative to the start codon (bp; aligned to Kp-M coordinates)", fontsize=7.5)

# two visual axes, both keyed: shape = differs or shared, fill = present or absent
H = [Line2D([], [], ls="none", marker="^", ms=5.4, mfc="#666", mec="#666", label="STRE"),
     Line2D([], [], ls="none", marker="D", ms=4.4, mfc="#666", mec="#666", label="YAP1"),
     Line2D([], [], ls="none", marker="|", ms=6.5, mew=1.1, color="#666", label="tick = shared"),
     Line2D([], [], ls="none", marker="^", ms=5.4, mfc="#666", mec="#666", label="filled = present in every copy"),
     Line2D([], [], ls="none", marker="^", ms=5.4, mfc="white", mec="#666", mew=1.1, label="open = absent"),
     Line2D([], [], color=MOUSE_C, lw=1.4, label="mouse, 7 strains"),
     Line2D([], [], color=GP_C, lw=1.4, label="guinea pig, 2 strains")]
_h = {h.get_label(): h for h in H}
_order = ["STRE", "tick = shared", "filled = present in every copy", "mouse, 7 strains",
          "YAP1", "", "open = absent", "guinea pig, 2 strains"]
_blank = Line2D([], [], ls="none", marker="none", label="")
ax.legend(handles=[_h.get(k, _blank) for k in _order], loc="upper left", bbox_to_anchor=(0.0, 1.0), ncol=4,
          frameon=False, fontsize=7.5, handlelength=1.5, columnspacing=1.6, handletextpad=0.5, borderaxespad=0)

for ext in ("svg", "png"):
    fig.savefig(f"{HERE}/stre_turnover.{ext}", dpi=600)
n_loss = sum(r["event"] == "loss" for r in rows); n_gain = sum(r["event"] == "gain" for r in rows)
assert n_loss + n_gain == N_DIFF
print(f"saved: {len(GENES)} genes, {n_loss} losses, {n_gain} gains, {N_SHARED} shared sites; "
      f"{drawn['symbol']} symbols and {drawn['tick']} ticks drawn")
print("  events per gene: " + ", ".join(f"{k} event{'s' if k > 1 else ''} in {v} gene{'s' if v > 1 else ''}"
                                        for k, v in sorted(DIST.items())) + f" = {N_DIFF}")
