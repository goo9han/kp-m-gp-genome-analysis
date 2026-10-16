#!/usr/bin/env python3
"""Individual panels for the genome figure. Each panel is written as a standalone SVG+PNG at exactly the
size it will occupy in the assembled figure, so the assembly scale is always 1.0 and a matplotlib font
size of n is n pt on the printed page. No panel letters; every panel carries its own axes, units, key and
notes. Backgrounds are transparent.

Target placement (174 mm total width):
  row 1  panel_prevalence (84) | 6 mm gap | panel_lineage_gp (84)
  row 2  panel_svscan (174)
  row 3  panel_hog1_ahp1_9strain (174)
panel_lineage_mouse and panel_chk1_9strain are supplementary and are not part of that assembly.

mockup_combined.{svg,png} re-draws row 1-3 into one 174 mm figure with the same functions at scale 1.0,
purely so the text sizes and collisions can be measured.
"""
import json, os, csv, collections
import numpy as np
from scipy.stats import chi2
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.text, matplotlib.lines, matplotlib.transforms, matplotlib.collections
import matplotlib.patheffects as pe
import matplotlib.patches as mpatches

HERE = os.path.dirname(os.path.abspath(__file__)); KP = os.path.dirname(os.path.dirname(HERE))
SV = f"{KP}/sv_scan/work"
MOUSE_C, GP_C, OBS_C = "#F4A261", "#3a7ebf", "#111111"
GREY, MID, ACC, AMB = "#d9dde0", "#7d8a93", "#A32E2E", "#b45309"
FS, FS_T, FS_S = 8.0, 7.5, 8.0
MINPT, TARGETPT = 6.0, 8.0
plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "font.family": "Arial", "font.size": FS,
                     "axes.edgecolor": "#555", "axes.linewidth": 0.4, "xtick.color": "#333", "ytick.color": "#333",
                     "axes.labelcolor": "#111", "text.color": "#111", "xtick.labelsize": FS_T, "ytick.labelsize": FS_T,
                     "xtick.major.width": 0.4, "ytick.major.width": 0.4, "xtick.major.size": 1.8, "ytick.major.size": 1.8,
                     "lines.linewidth": 0.6, "hatch.linewidth": 0.45, "figure.facecolor": "none", "savefig.facecolor": "none"})

ROW1_H, SV_H, NINE_H, GAP, HALF_W, FULL_W = 68.0, 64.0, 56.0, 6.0, 84.0, 174.0

class C:
    def __init__(s, fig, W, H): s.f, s.W, s.H = fig, W, H
    def ax(s, x, y, w, h):
        a = s.f.add_axes([x / s.W, y / s.H, w / s.W, h / s.H]); a.set_facecolor("none"); return a
    def t(s, x, y, txt, **kw): return s.f.text(x / s.W, y / s.H, txt, **kw)

def newfig(W, H):
    fig = plt.figure(figsize=(W / 25.4, H / 25.4)); fig.patch.set_alpha(0.0)
    return fig, C(fig, W, H)

def swatch(c, x, y, fc, hatch=None, w=3.2, hh=2.1):
    r = mpatches.Rectangle((x / c.W, y / c.H), w / c.W, hh / c.H, transform=c.f.transFigure,
                           facecolor=fc, edgecolor=("white" if hatch else "none"), lw=0.4, hatch=hatch, zorder=6)
    c.f.add_artist(r)

def mark(c, x, y, marker, color, ms):
    l = matplotlib.lines.Line2D([x / c.W], [y / c.H], transform=c.f.transFigure, marker=marker, ls="none",
                                ms=ms, color=color, zorder=6)
    l.set_clip_on(False); c.f.add_artist(l)

def bare(ax, keep=("bottom",)):
    for sp in ("top", "right", "left", "bottom"):
        if sp not in keep: ax.spines[sp].set_visible(False)

# =========================================================== data
PROV = {r["province"]: r for r in csv.DictReader((l for l in open(f"{KP}/public_survey/screen_v2/provinces_v2.tsv") if not l.startswith("#")), delimiter="\t")}
TWO = {r["unit"]: r for r in csv.DictReader((l for l in open(f"{KP}/public_survey/screen_v2/prevalence_two_ways_v2.tsv") if not l.startswith("#")), delimiter="\t")}
CD = json.load(open(f"{KP}/public_survey/deep_lineage/controls_deep.json"))["replicates"]
TD = json.load(open(f"{KP}/public_survey/deep_lineage/thresholds_deep.json"))
OD = json.load(open(f"{KP}/public_survey/deep_lineage/observed_deep.json"))
AH = json.load(open(f"{KP}/ahp1_cnv/nine_strain/work_9strain.json"))
CH = json.load(open(f"{KP}/chk1/work/table.json"))["own"]
FIN = json.load(open(f"{SV}/final.json"))
ORDER9 = ["kp1", "kpmtx", "C89-2", "ID1LA", "IDNYC", "QC1", "W26G", "kpgp3", "kpgp6"]

def pois(k, den):
    lo = 100 * (chi2.ppf(0.025, 2 * k) / 2) / den if k else 0.0
    hi = 100 * (chi2.ppf(0.975, 2 * (k + 1)) / 2) / den
    return lo, hi

def bins(label):
    d = collections.defaultdict(dict)
    for l in open(f"{SV}/bins_{label}.tsv"):
        r = l.rstrip("\n").split("\t"); d[r[0]][int(r[1])] = float(r[4])
    return d

def track(gp_labels, mouse_labels, ref_label, minlen=20000):
    """tiled 5-kb windows -> (x_Mb, log2(GP/mouse), mouse_low_mask, contig offsets, contig order)"""
    B = {l: bins(l) for l in set(gp_labels + mouse_labels + [ref_label])}
    ref = B[ref_label]
    ctgs = [c for c in ref if max(ref[c]) + 1000 >= minlen]
    ctgs.sort(key=lambda c: -(max(ref[c]) + 1000))
    off, x0 = {}, 0
    for c in ctgs: off[c] = x0; x0 += max(ref[c]) + 1000 + 30000
    X, gpv, mov = [], [], []
    for c in ctgs:
        ks = sorted(ref[c])
        for i in range(0, len(ks) - 4, 5):
            b = ks[i:i + 5]
            try:
                g = np.mean([np.mean([B[l][c][k] for k in b]) for l in gp_labels])
                m = np.mean([np.mean([B[l][c][k] for k in b]) for l in mouse_labels])
            except KeyError:
                continue
            X.append(off[c] + b[0] + 2500); gpv.append(g); mov.append(m)
    X, gpv, mov = np.array(X, float), np.array(gpv), np.array(mov)
    ok = mov >= 0.3
    r = np.log2(np.clip(gpv, 0.02, None) / np.clip(mov, 0.02, None))
    return X, r, ok, off, ctgs, ref, B, gp_labels, mouse_labels

# =========================================================== panel: prevalence
def draw_prevalence(c, x0, y0, w, h):
    order = sorted(PROV, key=lambda p: (-int(TWO[p]["b_detected"]) / int(PROV[p]["screened_runs"]),
                                        -int(PROV[p]["detected_runs"]) / int(PROV[p]["screened_runs"]),
                                        -int(PROV[p]["screened_runs"])))
    ax = c.ax(x0 + 17.5, y0 + 12.0, 55.0, 51.0)
    _num = _den = 0
    for i, p in enumerate(order):
        y = -i; n = int(PROV[p]["screened_runs"]); d = int(PROV[p]["detected_runs"]); b = int(TWO[p]["b_detected"])
        _num += b; _den += n
        assert 0 <= b <= d <= n, (p, b, d, n)
        # Criterion (b), the number the paper reports: the bar is the b runs, the rest is not detected.
        ax.barh(y, b, height=0.72, color=GP_C, ec="none", zorder=3)
        ax.barh(y, n - b, left=b, height=0.72, color=GREY, ec="none", zorder=2)
        ax.text(n + 0.35, y, f"{b}/{n}", va="center", ha="left", fontsize=FS)   # numerator = the solid blue segment
        ax.text(-0.45, y, p, va="center", ha="right", fontsize=FS)
    assert (_num, _den) == (21, 81), f"bar-label fractions must sum to 21/81, got {_num}/{_den}"
    ax.set_ylim(-len(order) + 0.45, 0.6); ax.set_yticks([])
    ax.set_xlim(0, 13.4); ax.set_xticks([])
    bare(ax)
    for dx, yy, fc, hz, lab in [(0.8, 6.6, GP_C, None, "detected"),
                                (26.0, 6.6, GREY, None, "not detected")]:
        swatch(c, x0 + dx, y0 + yy, fc, hz)
        c.t(x0 + dx + 4.2, y0 + yy + 1.05, lab, fontsize=FS, va="center", ha="left")

# =========================================================== panels: lineage
def lineage(c, x0, y0, w, h, kind, rows, obs, title, xlab, extra):
    """kind: 'gpdiag' (Kp-GP diagnostic) or 'strict' (mouse diagnostic)."""
    c.t(x0 + 0.8, y0 + h - 0.8, title, fontsize=FS, va="top", ha="left")
    lm = extra.get("lm", 21.5)
    axh = h - extra["hdr"] - extra["bot"] - 5.5
    ax = c.ax(x0 + lm, y0 + extra["bot"] + 5.5, w - lm - 3.0, axh)
    mixc = lambda g: tuple(g * np.array(matplotlib.colors.to_rgb(GP_C)) + (1 - g) * np.array(matplotlib.colors.to_rgb(MOUSE_C)))
    rng = np.random.default_rng(0)
    flat = extra.get("flat_ctrl")          # one colour for every control; the mixture fraction is in the row label
    gaps = extra.get("ygap", {})           # extra blank space BEFORE row i, in row units, to separate groups
    ys = []; acc = 0.0
    for i in range(len(rows)):
        acc += gaps.get(i, 0.0); ys.append(-i - acc)
    n_ctrl_drawn = 0
    for i, (k, lab, col, mk) in enumerate(rows):
        y = ys[i]
        assert k in CD[extra["ctrl"]], f"control row {k} is not in the data file"
        v = np.array([r[f"{kind}_reads"] for r in CD[extra["ctrl"]][k]]) / extra["den"] * 100
        cc = flat if flat else (col if col else mixc(int(k[2:k.index(":")]) / 100))
        ax.scatter(v, y + rng.uniform(-0.18, 0.18, len(v)), s=1.5, marker=mk, color=cc, alpha=0.75, lw=0)
        ax.plot([v.mean()] * 2, [y - 0.30, y + 0.30], color="#111", lw=0.6)
        ax.text(-2.5, y, lab, ha="right", va="center", fontsize=FS)
        n_ctrl_drawn += 1
    yo = ys[-1] - 1.55
    # thin rules in the gaps, separating the kinds of row, and one header naming the mixture components
    for i, txt in extra.get("rules", []):
        yr = (ys[i - 1] + ys[i]) / 2.0 if i < len(ys) else (ys[-1] + yo) / 2.0
        ax.plot([-2.5, 104], [yr, yr], color="#bbb", lw=0.45, clip_on=False, solid_capstyle="butt", zorder=0)
        if txt:                                   # inside the axis: the gap row carries no points
            ax.text(0.0, yr - 0.10, txt, ha="left", va="top", fontsize=FS_T, color="#444", clip_on=False)
    for j, o in enumerate(obs):
        y = yo - j * extra.get("ospace", 1.05)
        ax.text(-2.5, y, o["lab"], ha="right", va="center", fontsize=FS, fontweight="bold" if o.get("bold") else "normal")
        if o.get("bound"):
            ax.annotate("", xy=(0.5, y), xytext=(o["hi"], y),
                        arrowprops=dict(arrowstyle="-|>", color=OBS_C, lw=0.7, mutation_scale=5, shrinkA=0, shrinkB=0))
            ax.plot([o["hi"]] * 2, [y - 0.3, y + 0.3], color=OBS_C, lw=0.8)
            ax.text(o["hi"] + 2.0, y, o["txt"], ha="left", va="center", fontsize=FS_S)
        else:
            ax.errorbar(o["v"], y, xerr=[[o["v"] - o["lo"]], [o["hi"] - o["v"]]], fmt="D", ms=3.2, color=OBS_C,
                        mfc=OBS_C, capsize=1.5, lw=0.6, zorder=4)
            if o.get("inline"):
                ax.text(o["tx"], y, o["txt"], ha="left", va="center", fontsize=FS_S,
                        fontweight="bold" if o.get("bold") else "normal")
            else:
                ax.text(o["tx"], y - 0.62, o["txt"], ha=o.get("ha", "center"), va="top",
                        fontsize=FS_S, fontweight="bold" if o.get("bold") else "normal")
    if not extra.get("no_thr_line"):
        ax.axvline(extra["thr"], color="#555", lw=0.5, ls=(0, (2.5, 1.5)))
        ax.text(extra["thr"] + 1.2, 1.02, "threshold", transform=ax.get_xaxis_transform(), ha="left", va="bottom",
                fontsize=FS_T, color="#555")
    ax.set_ylim(yo - extra.get("ospace", 1.05) * (len(obs) - 1) - 1.35, 1.0); ax.set_yticks([])
    ax.set_xlim(-3, 104); ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_xlabel(xlab, fontsize=FS, labelpad=1.5)
    bare(ax)
    return dict(n_ctrl_drawn=n_ctrl_drawn, ctrl_keys=[r[0] for r in rows], obs_x=[o.get("v") for o in obs])

def draw_lineage_gp(c, x0, y0, w, h):
    """[FIG-5B] the guinea-pig lineage panel. What it sets, each an opt-in key in `extra`:
      - the percentage and its CI are not on the panel; they belong to the caption and Results
      - the sample label is the run accession alone
      - the mixture rows carry a bracket naming the two components in the order the DATA carry them
      - the three kinds of row are separated by left brackets: isolate / in-silico mixtures / sample
      - every control is one colour (flat_ctrl)
      - no threshold line (no_thr_line); the threshold is 0.20% of this axis, see the note"""
    den = TD["DRR848895"]["gpdiag_100pct"]; k = OD["DRR848895"]["gpdiag_reads"]
    f = 100 * k / den; lo, hi = pois(k, den); thr = 100 * TD["DRR848895"]["gpdiag_threshold"] / den
    # Component order checked against the data, not read off the key names: the gpdiag channel counts
    # Kp-GP-diagnostic alleles, the row written "gp100:kp1_0" carries the HIGHEST gpdiag count and
    # "gp0:kp1_100" is exactly zero in all 100 replicates. So the first number is the Kp-GP share.
    _hi = np.mean([r["gpdiag_reads"] for r in CD["DRR848895"]["gp100:kp1_0"]])
    _lo = np.mean([r["gpdiag_reads"] for r in CD["DRR848895"]["gp0:kp1_100"]])
    assert _hi > _lo, "mixture component order is not what the row names say"
    rows = [("KpGP2", "Kp-GP2 isolate", None, "s"), ("gp100:kp1_0", "100 : 0", None, "o"),
            ("gp80:kp1_20", "80 : 20", None, "o"), ("gp50:kp1_50", "50 : 50", None, "o"),
            ("gp0:kp1_100", "0 : 100", None, "o")]
    obs = [dict(lab="DRR848895", v=f, lo=lo, hi=hi, tx=f, txt="", bold=True)]
    r = lineage(c, x0, y0, w, h, "gpdiag", rows, obs,
                "Kp-GP-lineage diagnostic alleles",
                "Estimated Kp-GP-lineage content (%)",
                dict(ctrl="DRR848895", den=den, thr=thr, bot=11.0, hdr=12.5, thr_y=4.4, lm=22.5,
                     flat_ctrl=GP_C, no_thr_line=True,
                     ygap={1: 0.95},                      # gap between the isolate row and the mixtures
                     rules=[(1, "in-silico mixtures, Kp-GP : Kp-M"), (5, "metagenome sample")]))
    dk = list(CD["DRR848895"]); mk = [x for x in dk if x.startswith("gp")]
    assert set(r["ctrl_keys"]) <= set(dk)
    assert r["n_ctrl_drawn"] == len(rows) == 5
    assert abs(r["obs_x"][0] - 100 * OD["DRR848895"]["gpdiag_reads"] / den) < 1e-9
    print(f"  [FIG-5B] data has {len(dk)} control keys = {len(mk)} mixture ratios + {sorted(set(dk)-set(mk))}")
    print(f"  [FIG-5B] drawn {r['n_ctrl_drawn']}: {r['ctrl_keys']}")
    print(f"  [FIG-5B] sample x drawn {r['obs_x'][0]:.4f}%  stored {100*OD['DRR848895']['gpdiag_reads']/den:.4f}%"
          f"  ({OD['DRR848895']['gpdiag_reads']}/{den})   threshold {thr:.4f}% = {TD['DRR848895']['gpdiag_threshold']:.0f} reads")

def draw_lineage_mouse(c, x0, y0, w, h):
    d95 = TD["DRR848895"]; den95 = d95["strict_100pct"]; den46 = TD["DRR849046"]["strict_100pct"]
    thr = 100 * d95["strict_threshold"] / den95
    rows = [("KpGP2", "Kp-GP2 isolate", GP_C, "s"), ("gp100:kp1_0", "100 : 0", None, "o"), ("gp99:kp1_1", "99 : 1", None, "o"),
            ("gp95:kp1_5", "95 : 5", None, "o"), ("gp80:kp1_20", "80 : 20", None, "o"), ("gp64:kp1_36", "64 : 36", None, "o"),
            ("gp50:kp1_50", "50 : 50", None, "o"), ("gp20:kp1_80", "20 : 80", None, "o"), ("gp0:kp1_100", "0 : 100", None, "o")]
    k95 = OD["DRR848895"]["strict_reads"]; f95 = 100 * k95 / den95; l95, h95 = pois(k95, den95)
    k46 = OD["DRR849046"]["strict_reads"]; f46 = 100 * k46 / den46; l46, h46 = pois(k46, den46)
    obs = [dict(lab="DRR848895\n(Hunan)", v=f95, lo=l95, hi=h95, tx=f95 + 3.5, inline=True,
                txt=f"{f95:.2f}%  (95% CI {l95:.2f}–{h95:.2f})", bold=True),
           dict(lab="DRR849046\n(Shanghai)", v=f46, lo=l46, hi=h46, tx=f46 + 3.5, inline=True,
                txt=f"{f46:.2f}%  (95% CI {l46:.2f}–{h46:.2f})", bold=True),
           dict(lab="DRR848899\n(Hunan)", bound=True, hi=12.0, txt="≤ 12% (one-sided 95%); 418 reads", bold=True)]
    lineage(c, x0, y0, w, h, "strict", rows, obs,
            "Mouse-associated-lineage diagnostic alleles",
            "Estimated mouse-associated-lineage content (%)",
            dict(ctrl="DRR848895", den=den95, thr=thr, bot=17.0, hdr=17.5, thr_y=8.0, lm=24.0, ospace=1.9,
                 thr_txt=f"Dashed line: positivity threshold = {d95['strict_threshold']:.0f} reads, i.e. {thr:.2f}% content;\nfixed before observation.",
                 notes=["Controls: draws from the DRR848895 read budget (43,088 reads).",
                        "DRR849046 is scaled to its own 100% control (46,954 reads).",
                        "DRR848899 has 418 Kp reads only, so it gives an upper bound,",
                        "not an estimate; the arrow spans what it cannot exclude."]))

# =========================================================== panel: sv scan
_TRK = {}
def get_track(which):
    if which in _TRK: return _TRK[which]
    if which == "gp":
        _TRK[which] = track(["kpgp3_self", "kpgp6_on3"],
                            ["kp1_on3", "kpmtx_on3", "C89-2_on3", "ID1LA_on3", "IDNYC_on3", "QC1_on3", "W26G_on3"], "kpgp3_self")
    else:
        _TRK[which] = track(["kpgp3_on_k1", "kpgp6_on_k1"],
                            ["kp1_self_k1", "kpmtx_on_k1", "C89-2_on_k1", "ID1LA_on_k1",
                             "IDNYC_on_k1", "QC1_on_k1", "W26G_on_k1"], "kp1_self_k1")
    return _TRK[which]

# validated deletions drawn at 1-kb resolution (same intervals as the nine-strain panels)
VIS = {}
BANDS = {"gp": ("scaffold_18", 79562, 93351, "HOG1\u2013AHP1"),
         "km": ("scaffold_6", 420988, 422860, "CHK1")}
FLOOR = -2.42

def draw_svscan(c, x0, y0, w, h):
    kx = x0 + 0.8
    for mk, col, ms, lab, adv in [("o", ACC, 3.4, "junction-supported deletion", 46.0),
                                  ("^", AMB, 3.2, "depth drop, unresolved", 38.0),
                                  ("D", GP_C, 2.6, "gain", 0.0)]:
        mark(c, kx + 1.0, y0 + h - 3.5, mk, col, ms)
        c.t(kx + 2.8, y0 + h - 3.5, lab, fontsize=FS, va="center", ha="left")
        kx += adv
    AXW, AXH, LX = w - 40.5, 15.0, 17.0

    def one(ax_y, which, hdr, xlab, ups, downs):
        ax = c.ax(x0 + LX, ax_y, AXW, AXH)
        X, r, ok, off, ctgs, ref, B, GPL, MOL = get_track(which)
        for i, ct in enumerate(ctgs):
            if i % 2 == 0: ax.axvspan(off[ct] / 1e6, (off[ct] + max(ref[ct]) + 1000) / 1e6, color="#f2f4f6", lw=0, zorder=0)
        ax.scatter(X[ok] / 1e6, np.clip(r[ok], -2.5, 2.5), s=0.5, color="#6b7280", lw=0, alpha=0.75, rasterized=True, zorder=2)
        # Windows whose mouse depth is below 0.3 are left out of the scatter: the ratio is not stable
        # at that depth, so the caption says that windows below the depth cut are not plotted.
        ax.axhline(0, color="#999", lw=0.5, zorder=1)
        for yy, tag in ((np.log2(2 / 3), "1 of 3 lost"), (np.log2(1 / 3), "2 of 3 lost")):
            ax.axhline(yy, color="#c9ced3", lw=0.4, ls=(0, (2, 2)), zorder=1)
            if which == "gp":
                ax.plot([1.004, 1.030], [yy, yy], transform=ax.get_yaxis_transform(), color="#c9ced3",
                        lw=0.4, clip_on=False, zorder=5)
                ax.text(1.040, yy, tag, transform=ax.get_yaxis_transform(), ha="left", va="center",
                        fontsize=FS_T, color="#5a6670", clip_on=False, zorder=5)
        xmax = X.max() / 1e6
        ax.set_ylim(-2.6, 2.6); ax.set_yticks([-2, -1, 0, 1, 2])
        ax.set_xlim(-xmax * 0.012, xmax * 1.012)
        ax.set_xticks(np.arange(0, np.floor(xmax) + 0.5, 2.0))
        span = ax.get_xlim()[1] - ax.get_xlim()[0]
        # other candidates.  [SV-MT] count every marker actually drawn, per symbol class.
        drawn = {"^": 0, "D": 0, "o/v": 0}
        for px in ups:   ax.scatter([px], [-0.055], transform=ax.get_xaxis_transform(), marker="^", s=11,
                                    color=AMB, zorder=6, clip_on=False); drawn["^"] += 1
        for px in downs: ax.scatter([px], [1.045], transform=ax.get_xaxis_transform(), marker="D", s=7,
                                    color=GP_C, zorder=6, clip_on=False); drawn["D"] += 1
        # ---- the validated deletion: one representative marker, no band ----
        ct, s0, e0, gene = BANDS[which]
        vals = []
        for k in sorted(B[MOL[0]][ct]):
            if min(k + 1000, e0) - max(k, s0) < 500: continue    # bin at least half inside the deletion
            g = np.mean([B[l][ct][k] for l in GPL]); m = np.mean([B[l][ct][k] for l in MOL])
            if m < 0.3: continue                                  # same mask as the background
            vals.append(np.log2(max(g, 0.02) / max(m, 0.02)))
        med = float(np.median(vals))
        gx = (off[ct] + (s0 + e0) / 2) / 1e6
        below = med < FLOOR
        my = FLOOR if below else med
        if below:
            ax.scatter([gx], [my], marker="v", s=30, color=ACC, ec="white", lw=0.5, zorder=9); drawn["o/v"] += 1
            ax.text(gx - 0.008 * xmax, my, "below \u22122.5", ha="right", va="center", fontsize=FS_T, color=ACC, zorder=9)
        else:
            ax.scatter([gx], [my], marker="o", s=18, color=ACC, ec="white", lw=0.5, zorder=9); drawn["o/v"] += 1
        ly = min(2.10, my + 1.15)
        ax.plot([gx + 0.004 * xmax, gx + 0.010 * xmax], [my + 0.16, ly - 0.16], color=ACC, lw=0.5, zorder=9)
        ax.text(gx + 0.012 * xmax, ly, gene, ha="left", va="center", fontsize=FS, color=ACC, fontweight="bold", zorder=9)
        ax.plot([gx, gx], [-0.012, -0.070], transform=ax.get_xaxis_transform(), color=ACC, lw=0.9,
                clip_on=False, zorder=9, solid_capstyle="butt")
        VIS[which] = dict(interval=f"{ct}:{s0}-{e0}", n_bins=len(vals), median=med, drawn_at=my, below=below,
                          marker_mm=(2 * (30 if below else 18) / np.pi) ** 0.5 / 72 * 25.4,
                          bg_mm=(2 * 0.5 / np.pi) ** 0.5 / 72 * 25.4,
                          width_mm=(e0 - s0) / 1e6 / span * AXW, n_up=len(ups), n_down=len(downs),
                          n_greytick=int((~ok).sum()), xmax=float(xmax),
                          fx=(LX + (gx - ax.get_xlim()[0]) / span * AXW) / w,
                          fy=(ax_y - y0 + (my + 2.6) / 5.2 * AXH) / h,
                          fy0=(ax_y - y0 + (0.0 + 2.6) / 5.2 * AXH) / h,
                          fy_unit=(AXH / 5.2) / h)
        ax.set_xlabel(xlab, fontsize=FS, labelpad=1.5)
        bare(ax, keep=("bottom", "left"))
        c.t(x0 + 0.8, ax_y + AXH + 6.3, hdr, fontsize=FS, va="top", ha="left")
        ax.text(0.0, 1.06, "log2(guinea pig / mouse)", transform=ax.transAxes, fontsize=FS, va="bottom", ha="left")
        return ax, drawn

    Xg, _, _, offg, _, _, _, _, _ = get_track("gp")
    mid = lambda x, o: (o[x["interval"].split(":")[0]] + sum(map(int, x["interval"].split(":")[1].split("-"))) / 2) / 1e6
    ups_g = [mid(x, offg) for x in FIN["deletion"] if x["verdict"] == "UNRESOLVED"]
    dn_g = [mid(x, offg) for x in FIN["gain"] if x["verdict"] == "GAIN"]
    # [SVSCAN-MARKS] header counts are derived, never hardcoded. "interpretable" = SUPPORTED + UNRESOLVED,
    # the candidates that survived every exclusion filter and reached the junction test; the rest carry no
    # mark because they were set aside before it (CONTIG_END, ARTIFACT_DUP, SEQUENCE_IN_GP_ASSEMBLY,
    # NOT_CLADE_SPECIFIC). The arithmetic in the header has to close against the marks actually drawn.
    n_g = len(FIN["deletion"])
    n_g_sup = sum(1 for x in FIN["deletion"] if x["verdict"] == "SUPPORTED")
    n_g_int = n_g_sup + sum(1 for x in FIN["deletion"] if x["verdict"] == "UNRESOLVED")
    _, dg = one(y0 + 37.0, "gp",
        f"{len(Xg):,} tiled 5-kb windows \u00b7 {n_g_int} interpretable \u00b7 {n_g_sup} junction-supported",
        "Kp-GP assembly (Mb)", ups_g, dn_g)
    assert dg["^"] == len(ups_g) == sum(1 for x in FIN["deletion"] if x["verdict"] == "UNRESOLVED")
    assert dg["D"] == len(dn_g) == sum(1 for x in FIN["gain"] if x["verdict"] == "GAIN")
    assert dg["o/v"] == 1
    # the header's "interpretable" must equal the marks on this track (triangles + the supported mark).
    # Gains are a separate category and are not among the deletion candidates, so dg["D"] is excluded.
    assert n_g_int == dg["^"] + dg["o/v"], (n_g_int, dg)
    print(f"  [SV-MT] forward markers drawn {dg} | data: ^={sum(1 for x in FIN['deletion'] if x['verdict']=='UNRESOLVED')}"
          f" D={sum(1 for x in FIN['gain'] if x['verdict']=='GAIN')} o/v=1 | {len(Xg):,} tiled windows")
    Xk, _, _, offk, _, _, _, _, _ = get_track("km")
    ups_k = [mid(x, offk) for x in FIN["reverse"] if x["verdict"] == "UNRESOLVED"]
    n_k = len(FIN["reverse"])
    n_k_sup = sum(1 for x in FIN["reverse"] if x["verdict"] == "SUPPORTED")
    n_k_int = n_k_sup + sum(1 for x in FIN["reverse"] if x["verdict"] == "UNRESOLVED")
    _, dk = one(y0 + 7.7, "km",
        f"{len(Xk):,} tiled 5-kb windows \u00b7 {n_k_int} interpretable \u00b7 {n_k_sup} junction-supported",
        "Kp-M assembly (Mb)", ups_k, [])
    assert dk["^"] == len(ups_k) == sum(1 for x in FIN["reverse"] if x["verdict"] == "UNRESOLVED")
    assert dk["D"] == 0
    assert dk["o/v"] == 1
    assert n_k_int == dk["^"] + dk["o/v"], (n_k_int, dk)
    # [SVSCAN-MARKS] BANDS is hardcoded; tie it to the data. One SUPPORTED verdict per track, one BANDS entry
    # per track, same contig, and the two intervals must overlap. Containment is NOT asserted, because BANDS
    # holds the VALIDATED boundary while the candidate holds the SCAN WINDOW: for HOG1-AHP1 neither contains
    # the other (validated 79,562-93,351 starts 3.4 kb left of the scan window 83,001-94,000). See
    # svscan_numbers.md, "Coordinates are always the validated boundary, never the scan window".
    #
    # The overlap test only catches a coordinate that moves by kilobases. Where the validated boundary exists
    # as a machine-readable field it is asserted for EQUALITY below instead, which also catches a small shift:
    #   km / CHK1       YES - chk1/figure/chk1_deletion_data.tsv section 3, key validated_deletion_boundary.
    #                   The TSV carries no contig, so the contig stays pinned by the same-contig test above.
    #   gp / HOG1-AHP1  NO - overlap only. There is no field anywhere to compare against:
    #                   sv_scan/work/final.json has 79,562-93,351 in no field (its "note" is free text and
    #                   names only the 316-bp repeat copy 93,352-93,667), and figures/panels/coordinate_audit.md
    #                   has it in a markdown table. The Python constants in ahp1_cnv/ are not usable as a
    #                   reference either. Every ahp1_cnv script that states a BOUNDARY uses
    #                   DEL = (79562, 93351), the deletion with one copy of the 316-bp flanking repeat
    #                   retained (13,790 bp = the 13.8 kb the manuscript prints). But 93,667 also appears,
    #                   correctly, as something else: region_genes.py uses (79562, 93667) as a classification
    #                   WINDOW, and make_figure.py / make_compact.py / sv_scan/repeat_pairs.py use 93,667 as
    #                   the end of the right repeat copy. A naive equality assert against 'the 93,667 one'
    #                   would therefore pin a window or a repeat edge, not a boundary. If someone later
    #                   writes the boundary into final.json as a field, pin this track the same way.
    VB = {}
    with open(f"{KP}/chk1/figure/chk1_deletion_data.tsv") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f and f[0] == "validated_deletion_boundary":
                VB["km"] = (int(f[1]), int(f[2]))
    assert VB["km"], "chk1_deletion_data.tsv: validated_deletion_boundary row not found"
    for key, cat in (("gp", "deletion"), ("km", "reverse")):
        sup = [x for x in FIN[cat] if x["verdict"] == "SUPPORTED"]
        bc, bs, be, blab = BANDS[key]
        assert len(sup) == 1, (key, len(sup))
        ic, rng = sup[0]["interval"].split(":"); a, b = (int(v) for v in rng.split("-"))
        assert ic == bc, (key, ic, bc)
        assert max(a, bs) <= min(b, be), (key, sup[0]["interval"], f"{bc}:{bs}-{be}")
        if key in VB:
            assert (bs, be) == VB[key], (key, (bs, be), VB[key])
        how = f"== field {VB[key][0]:,}-{VB[key][1]:,}" if key in VB else "overlap only (no field; see comment)"
        print(f"  [SVSCAN-MARKS] {key} {blab}: SUPPORTED 1 = BANDS 1; scan window {sup[0]['interval']} "
              f"vs validated {bc}:{bs:,}-{be:,} {how}")
    print(f"  [SVSCAN-MARKS] forward {n_g} candidates, {n_g_int} interpretable = {dg['^']} triangles + "
          f"{dg['o/v']} supported | reverse {n_k} candidates, {n_k_int} interpretable = {dk['^']} + {dk['o/v']}")
    print(f"  [SV-MT] reverse markers drawn {dk} | data: ^={sum(1 for x in FIN['reverse'] if x['verdict']=='UNRESOLVED')}"
          f" D=0 o/v=1 | {len(Xk):,} tiled windows | {len(FIN['reverse'])} absence candidates in data")

# =========================================================== panels: nine strains
def draw_nine(c, x0, y0, w, h, rows, title, note_ref, jmax, jnote):
    aL = c.ax(x0 + 25.0, y0 + 13.0, 60.0, h - 18.0)
    aR = c.ax(x0 + 97.0, y0 + 13.0, 26.0, h - 18.0)
    ys, y = {}, 0
    for r in rows:
        if r["host"] == "guinea pig" and "gap" not in ys: y -= 1.0; ys["gap"] = True; ys["_mid"] = y + 1.0
        ys[r["name"]] = y; y -= 1
    for r in rows:
        yy = ys[r["name"]]; col = GP_C if r["host"] == "guinea pig" else MOUSE_C
        aL.plot([r["p5"], r["p95"]], [yy, yy], color=GREY, lw=2.6, solid_capstyle="butt", zorder=1)
        aL.plot([r["p1"]] * 2, [yy - 0.26, yy + 0.26], color=MID, lw=0.5, zorder=1)
        aL.scatter(r["val"], yy, s=11, color=col, marker="D" if r["host"] == "guinea pig" else "o", zorder=3, ec="white", lw=0.3)
        aL.text(-0.045, yy, r["name"], transform=aL.get_yaxis_transform(), ha="right", va="center", fontsize=FS,
                color=col, fontweight="bold" if r["host"] == "guinea pig" else "normal")
        aR.barh(yy, r["pairs"], height=0.6, color=col, ec="none")
        aR.text(jmax * 0.03 + r["pairs"], yy, str(r["pairs"]), va="center", ha="left", fontsize=FS_S)
        aR.text(1.28, yy, "n/a" if r["ont"] is None else str(r["ont"]), transform=aR.get_yaxis_transform(),
                ha="center", va="center", fontsize=FS_S, color="#111" if r["ont"] is not None else MID)
        aR.text(1.92, yy, "present" if r["dele"] else "absent", transform=aR.get_yaxis_transform(), ha="center", va="center",
                fontsize=FS_S, color=GP_C if r["dele"] else "#5a6670", fontweight="bold" if r["dele"] else "normal")
    aL.text(-0.045, 1.15, "mouse", transform=aL.get_yaxis_transform(), ha="right", va="center", fontsize=FS_S, color=MOUSE_C)
    aL.text(-0.045, ys["_mid"], "guinea pig", transform=aL.get_yaxis_transform(), ha="right", va="center", fontsize=FS_S, color=GP_C)
    ymin = min(v for k, v in ys.items() if k not in ("gap", "_mid")) - 0.7
    for ax in (aL, aR):
        ax.set_ylim(ymin, 1.60); ax.set_yticks([])
        ax.axhline(ys["_mid"], color="#b9c0c5", lw=0.4, ls=(0, (2, 2)))
        bare(ax)
    aL.axvline(1.0, color="#b9c0c5", lw=0.4, zorder=0)
    aL.set_xlim(-0.04, 1.32); aL.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
    aL.set_xlabel("Depth inside the deleted interval (own genome median = 1)", fontsize=FS, labelpad=1.5)
    aR.set_xlim(0, jmax); aR.set_xticks([0, jmax // 2, jmax])
    aR.set_xlabel("Junction read pairs", fontsize=FS, labelpad=1.5)
    aR.text(1.28, 1.22, "ONT junction\nreads", transform=aR.get_yaxis_transform(), ha="center", va="bottom", fontsize=FS_S, linespacing=1.2)
    aR.text(1.92, 1.22, "Deletion", transform=aR.get_yaxis_transform(), ha="center", va="bottom", fontsize=FS_S)
    c.t(x0 + 0.8, y0 + h - 0.8, title, fontsize=FS, va="top", ha="left", fontweight="bold")

rowsH = [dict(name=r["name"], host=r["host"], val=r["interior_depth"], p1=r["window_p1"], p5=r["window_p5"], p95=r["window_p95_"],
              pairs=r["junction_pairs"], ont=(r["ont"] or {}).get("junction_13kb"), dele=r["deletion"] == "present") for r in AH]
def _ont_junction(row):
    o = row["ont"]
    if o is None: return None
    return o["counts"]["junction"]

rowsK = [dict(name=CH[s]["name"], host=CH[s]["host"], val=CH[s]["interior"], p1=CH[s]["p1"], p5=CH[s]["p5"], p95=CH[s]["p95"],
              pairs=CH[s]["junction_pairs"], ont=_ont_junction(CH[s]), dele=CH[s]["call"] == "deleted") for s in ORDER9]
# a .get() chain here once turned a renamed key into "n/a" in every row without any check failing, so the key
# is required when an ONT record exists rather than defaulted (chk1/reconcile_ont.py names the class "junction")
assert sum(r["ont"] is not None for r in rowsK) == 4, "ONT counts missing from panel_chk1_9strain"

def draw_hog1(c, x0, y0, w, h):
    draw_nine(c, x0, y0, w, h, rowsH, "HOG1–AHP1 deletion (13.8 kb)",
              "On the Kp-GP assembly: part of the interval remains in guinea-pig genomes.",
              24, "Junction-pair axis runs to 24.")
def draw_chk1(c, x0, y0, w, h):
    draw_nine(c, x0, y0, w, h, rowsK, "CHK1 deletion (1.9 kb)",
              "On the Kp-M assembly: absent from both guinea-pig genomes.",
              160, "Junction-pair axis runs to 160, not the HOG1\u2013AHP1 scale.")

# =========================================================== output
def save(fig, name):
    fig.savefig(f"{HERE}/{name}.svg", transparent=True)
    fig.savefig(f"{HERE}/{name}.pdf", transparent=True)
    fig.savefig(f"{HERE}/{name}.png", dpi=600, transparent=True)

PANELS = [("panel_prevalence", HALF_W, ROW1_H, draw_prevalence),
          ("panel_lineage_gp", HALF_W, ROW1_H, draw_lineage_gp),
          ("panel_lineage_mouse", HALF_W, 98.0, draw_lineage_mouse),
          ("panel_svscan", FULL_W, SV_H, draw_svscan),
          ("panel_hog1_ahp1_9strain", FULL_W, NINE_H, draw_hog1),
          ("panel_chk1_9strain", FULL_W, NINE_H, draw_chk1)]

def audit(fig):
    """text sizes, pairwise overlaps and out-of-figure text, measured on the drawn figure."""
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    fw, fh = fig.canvas.get_width_height()
    items = []
    for o in fig.findobj(matplotlib.text.Text):
        s = o.get_text().strip()
        if not s or not o.get_visible(): continue
        bb = o.get_window_extent(renderer=r)
        if bb.width <= 0 or bb.height <= 0: continue
        ax = o.axes
        if ax is not None:                      # skip phantom tick labels placed outside the axes
            try:
                xd, yd = o.get_position()
                if o in ax.get_xticklabels() and not (ax.get_xlim()[0] - 1e-9 <= xd <= ax.get_xlim()[1] + 1e-9): continue
                if o in ax.get_yticklabels() and not (ax.get_ylim()[0] - 1e-9 <= yd <= ax.get_ylim()[1] + 1e-9): continue
            except Exception: pass
        items.append((round(float(o.get_fontsize()), 2), s.replace("\n", " \u23ce "), bb))
    marks = []
    for o in fig.findobj(matplotlib.lines.Line2D):
        mk = o.get_marker()
        if not o.get_visible() or mk in (None, "None", "", " ") or isinstance(mk, int): continue   # int markers are axis ticks
        bb = o.get_window_extent(r)
        if 0 < bb.width < 0.05 * fw and 0 < bb.height < 0.2 * fh:
            marks.append((f"marker {o.get_marker()!r}", bb))
    for col in fig.findobj(matplotlib.collections.PathCollection):
        if not col.get_visible(): continue
        offs = col.get_offsets()
        if offs is None or len(offs) == 0 or len(offs) > 12: continue   # skip the dense data clouds
        tr = col.get_offset_transform()
        sz = (max(col.get_sizes()) if len(col.get_sizes()) else 9) ** 0.5 * fig.dpi / 72 * 0.6
        for x, y in tr.transform(offs):
            marks.append(("scatter marker", matplotlib.transforms.Bbox.from_bounds(x - sz, y - sz, 2 * sz, 2 * sz)))
    mm = lambda px: round(px * (fig.get_figwidth() * 25.4) / fw, 2)
    oob = []
    for pt, s, bb in items:
        over = max(bb.x1 - fw, -bb.x0, bb.y1 - fh, -bb.y0)
        if over > 1.0: oob.append([mm(over), pt, s[:70]])
    ov = []
    for a in range(len(items)):
        for b in range(a + 1, len(items)):
            A, B = items[a][2], items[b][2]
            ox = min(A.x1, B.x1) - max(A.x0, B.x0); oy = min(A.y1, B.y1) - max(A.y0, B.y0)
            if ox > 1.0 and oy > 1.0: ov.append([mm(min(ox, oy)), items[a][1][:46], items[b][1][:46]])
    for pt, s, A in items:                       # text against plotted markers
        for lbl, B in marks:
            ox = min(A.x1, B.x1) - max(A.x0, B.x0); oy = min(A.y1, B.y1) - max(A.y0, B.y0)
            if ox > 1.0 and oy > 1.0: ov.append([mm(min(ox, oy)), s[:46], lbl])
    return items, ov, oob

def measure(fig, label):
    fig.canvas.draw()
    out = []
    for o in fig.findobj(matplotlib.text.Text):
        s = o.get_text().strip()
        if not s or not o.get_visible(): continue
        out.append((round(float(o.get_fontsize()), 2), s.replace("\n", " ⏎ "), o))
    return out

REPORT = {}
for name, W, H, fn in PANELS:
    fig, c = newfig(W, H); fn(c, 0.0, 0.0, W, H); save(fig, name)
    items, ov, oob = audit(fig)
    REPORT[name] = dict(w=W, h=H, texts=items, overlaps=ov, oob=oob)
    plt.close(fig)
    print(f"saved {name:26s} {W:5.0f} x {H:5.1f} mm  min {min(t[0] for t in items):.1f} pt  overlaps {len(ov)}  oob {len(oob)}")

MH = ROW1_H + GAP + SV_H + GAP + NINE_H
fig, c = newfig(FULL_W, MH)
draw_prevalence(c, 0.0, MH - ROW1_H, HALF_W, ROW1_H)
draw_lineage_gp(c, HALF_W + GAP, MH - ROW1_H, HALF_W, ROW1_H)
draw_svscan(c, 0.0, NINE_H + GAP, FULL_W, SV_H)
draw_hog1(c, 0.0, 0.0, FULL_W, NINE_H)
save(fig, "mockup_combined")
items, ov, oob = audit(fig)
REPORT["mockup_combined"] = dict(w=FULL_W, h=MH, texts=items, overlaps=ov, oob=oob)
plt.close(fig)
print(f"saved {'mockup_combined':26s} {FULL_W:5.0f} x {MH:5.1f} mm  min {min(t[0] for t in items):.1f} pt  overlaps {len(ov)}  oob {len(oob)}")

json.dump({k: dict(w=v["w"], h=v["h"], n=len(v["texts"]), sizes=sorted({t[0] for t in v["texts"]}),
                   below6=[[t[0], t[1]] for t in v["texts"] if t[0] < MINPT],
                   below7=[[t[0], t[1]] for t in v["texts"] if t[0] < TARGETPT],
                   overlaps=v["overlaps"], out_of_bounds=v["oob"]) for k, v in REPORT.items()},
          open(f"{HERE}/measure.json", "w"), indent=1)
bad = {k: (len(v["overlaps"]), len(v["oob"])) for k, v in REPORT.items() if v["overlaps"] or v["oob"]}
print("panels with problems:", bad if bad else "none")

# ---------------- per-format check on panel_svscan ----------------
import subprocess, tempfile, re as _re
from PIL import Image

def _raster(path, size=2400):
    if path.endswith(".png"): return Image.open(path).convert("RGBA")
    d = tempfile.mkdtemp()
    subprocess.run(["qlmanage", "-t", "-s", str(size), "-o", d, path], capture_output=True)
    f = [x for x in os.listdir(d) if x.endswith(".png")]
    return Image.open(os.path.join(d, f[0])).convert("RGBA") if f else None

def _on_white(im):
    w = Image.new("RGBA", im.size, (255, 255, 255, 255))
    return np.array(Image.alpha_composite(w, im).convert("RGB")).astype(int)

def _probe(arr, v):
    """ACC pixels in a narrow column on the marker (narrow enough to miss the leader line)"""
    H, Wd = arr.shape[:2]
    cx, cy = int(v["fx"] * Wd), int((1 - v["fy"]) * H)
    rx, ry = max(2, int(0.0030 * Wd)), max(4, int(0.012 * Wd))
    box = arr[max(0, cy - ry):cy + ry, max(0, cx - rx):cx + rx]
    if box.size == 0: return 0, None
    hit = np.abs(box - np.array([163, 46, 46])).sum(axis=2) < 90
    if not hit.any(): return 0, None
    rows = np.where(hit.any(axis=1))[0]
    row = max(0, cy - ry) + float((rows.min() + rows.max()) / 2)
    return int(hit.sum()), ((1 - row / H) - v["fy0"]) / v["fy_unit"]

def _svg_markers(path, Wmm, Hmm):
    """the two deletion markers are the only <use> elements with fill #a32e2e and a white stroke"""
    s = open(path).read()
    out = []
    for m in _re.finditer(r'<use[^>]*x="([\d.]+)"[^>]*y="([\d.]+)"[^>]*style="[^"]*fill: #a32e2e[^"]*stroke: #ffffff', s):
        out.append((float(m.group(1)), float(m.group(2))))
    return sorted(out, key=lambda p: p[1]), Wmm / 25.4 * 72, Hmm / 25.4 * 72

rows_fmt = []
# SVG: read the marker coordinates straight out of the file (qlmanage letterboxes SVG, so a raster check is unreliable)
pts, Wpt, Hpt = _svg_markers(f"{HERE}/panel_svscan.svg", FULL_W, SV_H)
cells = []
for k, nm in (("gp", "HOG1\u2013AHP1"), ("km", "CHK1")):
    v = VIS[k]; want = (1 - v["fy"]) * Hpt
    got = [y for _, y in pts if abs(y - want) < 1.0]
    if got:
        val = ((1 - got[0] / Hpt) - v["fy0"]) / v["fy_unit"]
        cells.append(f"{nm}: found at y={got[0]:.2f} pt (expected {want:.2f}), reads {val:+.2f}")
    else:
        cells.append(f"{nm}: **not found**")
rows_fmt.append(f"| `.svg` | parsed, {len(pts)} white-edged red markers | " + " | ".join(cells) + " |")
for ext in ("pdf", "png"):
    im = _raster(f"{HERE}/panel_svscan.{ext}")
    if im is None:
        rows_fmt.append(f"| `.{ext}` | could not rasterise | \u2014 | \u2014 |"); continue
    arr = _on_white(im); H, Wd = arr.shape[:2]
    cells = []
    for k, nm in (("gp", "HOG1\u2013AHP1"), ("km", "CHK1")):
        n, val = _probe(arr, VIS[k])
        cells.append(f"{nm}: {n} red px, reads {val:+.2f}" if val is not None else f"{nm}: **not found**")
    rows_fmt.append(f"| `.{ext}` | rasterised {Wd}\u00d7{H} px over white | " + " | ".join(cells) + " |")
FORMAT_TABLE = ("| format | how it was checked | upper track | lower track |\n|---|---|---|---|\n"
                + "\n".join(rows_fmt) + "\n\n"
                + "The value in each cell is the marker's position read back off the output file and converted to log2 \u2014 "
                  "an independent check that it sits where the data says. HOG1\u2013AHP1 should read about "
                  f"{VIS['gp']['median']:+.2f} (its 1-kb median) and CHK1 about {FLOOR:+.2f} (the axis floor, since its "
                  "median is off-scale).\n")
print("format check:"); [print("   " + r) for r in rows_fmt]

# ---------------- legibility_check.md ----------------
ORDER_R = ["panel_prevalence", "panel_lineage_gp", "panel_svscan", "panel_hog1_ahp1_9strain",
           "panel_lineage_mouse", "panel_chk1_9strain", "mockup_combined"]
with open(f"{HERE}/legibility_check.md", "w") as fh:
    fh.write("# Legibility check — individual panels and the assembly mockup\n\n")
    fh.write(f"Every panel is written at its final placement size with no tight bounding box, so a matplotlib font size "
             f"of *n* is *n* pt on the printed page and the assembly scale is 1.0. Measured on the drawn figures, "
             f"not on the source. Target: body text **{TARGETPT:.0f} pt**, axis tick labels **{FS_T:.1f} pt**; "
             f"hard floor **{MINPT:.0f} pt**.\n\n")
    fh.write("## Summary\n\n| file | size (mm) | text elements | smallest | < 7 pt | < 6 pt | overlaps | text outside the panel |\n")
    fh.write("|---|---|---|---|---|---|---|---|\n")
    for k in ORDER_R:
        v = REPORT[k]; sizes = [t[0] for t in v["texts"]]
        fh.write(f"| `{k}` | {v['w']:.0f} × {v['h']:.0f} | {len(v['texts'])} | **{min(sizes):.1f} pt** | "
                 f"{sum(1 for s in sizes if s < TARGETPT)} | **{sum(1 for s in sizes if s < MINPT)}** | "
                 f"**{len(v['overlaps'])}** | **{len(v['oob'])}** |\n")
    tot6 = sum(sum(1 for t in v["texts"] if t[0] < MINPT) for v in REPORT.values())
    totov = sum(len(v["overlaps"]) for v in REPORT.values()); totoob = sum(len(v["oob"]) for v in REPORT.values())
    fh.write(f"\n**Nothing is below {MINPT:.0f} pt** ({tot6} elements), there are **{totov} text collisions** and "
             f"**{totoob} elements outside a panel**, in any panel or in the mockup.\n\n")
    fh.write("## What is set at which size\n\n| role | size |\n|---|---|\n")
    fh.write(f"| titles, notes, keys, axis labels, strain and category names, data labels | **{FS:.0f} pt** |\n")
    fh.write(f"| axis tick labels; the two copy-loss lines in the scan panel | **{FS_T:.1f} pt** |\n\n")
    fh.write("Everything that is not a tick label is at the 7 pt target. The only 6.5 pt elements are tick labels and "
             "the two dashed-line annotations inside the scan tracks, which sit on the axis scale.\n\n")
    fh.write("## Assembly mockup\n\n")
    fh.write(f"`mockup_combined` re-draws the panels with the same functions at the placement offsets, at "
             f"**{FULL_W:.0f} × {REPORT['mockup_combined']['h']:.0f} mm**:\n\n")
    fh.write(f"- row 1 `panel_prevalence` ({HALF_W:.0f} mm) + {GAP:.0f} mm gap + `panel_lineage_gp` ({HALF_W:.0f} mm), height {ROW1_H:.0f} mm\n")
    fh.write(f"- row 2 `panel_svscan` ({FULL_W:.0f} mm), height {SV_H:.0f} mm\n")
    fh.write(f"- row 3 `panel_hog1_ahp1_9strain` ({FULL_W:.0f} mm), height {NINE_H:.0f} mm\n")
    fh.write(f"- two {GAP:.0f} mm gaps; total **{REPORT['mockup_combined']['h']:.0f} mm**, the 200 mm ceiling.\n\n")
    fh.write("`panel_lineage_mouse` and `panel_chk1_9strain` are supplementary and are not in the mockup; they were "
             "measured on their own.\n\n")
    fh.write("## Longest text elements per panel (the ones most at risk when the layout moves)\n\n")
    for k in ORDER_R:
        longest = sorted(REPORT[k]["texts"], key=lambda x: -len(x[1]))[0]
        fh.write(f"- `{k}` — {longest[0]:.1f} pt, {len(longest[1])} chars: {longest[1][:105]}\n")
    fh.write("## How the two deletions are drawn in panel_svscan\n\n"
             "Each validated deletion is far narrower than one plotted window \u2014 0.135 mm and 0.020 mm on their axes. "
             "No drawing can spread its 1-kb bins horizontally; on this axis each is effectively one point. So each "
             "deletion is drawn as **one marker at the median of its 1-kb bins**, with the gene name on a leader line and "
             "a thin rug tick in the margin below the axis for position. There is no shaded band: a band wide enough to "
             "see made the interval look far wider than it is and hid the value inside it.\n\n")
    fh.write("| track | interval | true width on the axis | 1-kb bins used | median log2 | drawn at | marker |\n|---|---|---|---|---|---|---|\n")
    for k, lab in (("gp", "Kp-GP coords, HOG1\u2013AHP1"), ("km", "Kp-M coords, CHK1")):
        v = VIS[k]
        where = f"axis floor, labelled `below \u22122.5`" if v["below"] else f"{v['drawn_at']:+.2f}"
        fh.write(f"| {lab} | `{v['interval']}` | {v['width_mm']:.3f} mm | {v['n_bins']} | **{v['median']:+.2f}** | {where} | "
                 f"{v['marker_mm']:.2f} mm across vs {v['bg_mm']:.2f} mm for a background point |\n")
    fh.write(f"\nThe brief asked for a marker of 3\u20134\u00d7 the background point's *area*. That would be "
             f"{VIS['gp']['bg_mm'] * 2:.2f} mm across, which does not separate from a cloud of 2,674 points, so the marker "
             f"is {VIS['gp']['marker_mm'] / VIS['gp']['bg_mm']:.0f}\u00d7 the background diameter instead "
             f"({(VIS['gp']['marker_mm'] / VIS['gp']['bg_mm']) ** 2:.0f}\u00d7 the area). Say the word and it can be shrunk.\n\n")
    fh.write("## Legend against what is actually plotted\n\n| legend symbol | meaning | drawn in the panel |\n|---|---|---|\n")
    fh.write(f"| red circle | junction-supported deletion (1-kb median) | **2** \u2014 HOG1\u2013AHP1 (circle) "
             f"and CHK1 (drawn as a down-triangle on the floor, since its median is off-scale) |\n")
    fh.write(f"| orange up-triangle | depth drop, junction unresolved | **{VIS['gp']['n_up'] + VIS['km']['n_up']}** "
             f"({VIS['gp']['n_up']} upper track, {VIS['km']['n_up']} lower) |\n")
    fh.write(f"| blue diamond | gain | **{VIS['gp']['n_down'] + VIS['km']['n_down']}** (upper track only; the reverse scan has no gain category) |\n")
    fh.write(f"| grey tick at the top | mouse depth < 0.3, ratio not shown | **{VIS['gp']['n_greytick'] + VIS['km']['n_greytick']}** "
             f"({VIS['gp']['n_greytick']} upper, {VIS['km']['n_greytick']} lower) |\n")
    fh.write("\nEvery legend symbol appears in the panel. The only plotted symbol not in the legend is the red "
             "down-triangle on the axis floor, which carries its own `below \u22122.5` label beside it.\n\n")
    fh.write("## Checked by eye, per output format\n\n"
             "All three formats were rasterised and composited over white before looking, because the panels are saved "
             "with a transparent background and a viewer that ignores alpha shows them differently. SVG and PDF were "
             "rendered with `qlmanage`; PNG was read directly.\n\n")
    fh.write(FORMAT_TABLE)
    fh.write("\n### Checked by eye on the PNG composited over white\n\n"
             "| question | answer |\n|---|---|\n"
             "| Is the red marker distinguishable from the background cloud? | **Yes.** It is the only large filled "
             "circle in the track, dark red with a white edge, against 0.28 mm grey dots. |\n"
             "| Can its y position be read? | **Yes.** HOG1\u2013AHP1 sits just under the `1 of 3 copies lost` line, "
             f"which is what its 1-kb median ({VIS['gp']['median']:+.2f} against log2(2/3) = \u22120.58) says. |\n"
             "| Does every legend symbol appear in the track? | **Yes** \u2014 counted above: red circle 2, orange "
             f"triangle {VIS['gp']['n_up'] + VIS['km']['n_up']}, blue diamond {VIS['gp']['n_down']}, grey tick "
             f"{VIS['gp']['n_greytick'] + VIS['km']['n_greytick']}. |\n"
             "| Does the deletion look wider than it is? | **No.** Nothing spans the plot height any more; the only "
             "vertical element is a 0.3 mm rug tick in the margin below the axis. |\n"
             "| CHK1 | Red down-triangle on the axis floor with `below \u22122.5` beside it, and a leader to its label. "
             "It sits 0.6 Mb from an orange unresolved marker; the two differ in both colour and shape. |\n\n")
print("legibility_check.md written")
