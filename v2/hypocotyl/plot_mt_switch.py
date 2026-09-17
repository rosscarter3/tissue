#!/usr/bin/env python3
"""The microtubule switch: simulation against the measurements.

Reads a print-flag-0 run of the hook model with InnerTissue::GrowingCore and
produces figures/fig7_mt_switch.png.

The predicted CMT axis is the maximal principal stress direction computed by
the TRBS solver (cell variables 0-2, with the stress anisotropy on variable 3).
Panel B uses the paper's own signed convention from hook_part.py,

    signed anisotropy = +a if hoop dominates, -a if longitudinal dominates

so the switch is a zero crossing.

Usage: python3 plot_mt_switch.py <run.out> [t_end_hours]
"""
import math, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import analyze as AB
import expdata as ED

OUT = "figures"
os.makedirs(OUT, exist_ok=True)
C_IN, C_OUT, C_EXP, C_SIM = "#15706c", "#b04d22", "#22272e", "#46337e"
plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.grid": True,
    "grid.alpha": 0.16, "grid.linewidth": 0.6, "legend.frameon": False,
})

NC = AB.NCIRC
INNER, OUTER = [0, NC - 1], [NC // 2 - 1, NC // 2]


def series(path, t_end):
    fr = AB.parse(path)
    n = len(fr)
    _, n_print = AB.print_times(t_end)
    t = np.array([t_end * i / n_print for i in range(n)])
    out = {k: [] for k in ("angle", "mt_in", "mt_out", "a_in", "a_out")}
    for v, c in fr:
        mi, ai, _, _ = AB.flank_mt(v, c, INNER)
        mo, ao, _, _ = AB.flank_mt(v, c, OUTER)
        out["angle"].append(AB.hook_angle(v))
        out["mt_in"].append(mi);  out["mt_out"].append(mo)
        # paper's sign convention: + circumferential, - longitudinal
        out["a_in"].append(ai if mi > 45 else -ai)
        out["a_out"].append(ao if mo > 45 else -ao)
    return t, {k: np.array(v) for k, v in out.items()}, fr


def crossing(t, mt):
    """First time the flank passes 45 deg, linearly interpolated."""
    for i in range(1, len(mt)):
        if mt[i - 1] > 45 >= mt[i]:
            f = (mt[i - 1] - 45) / (mt[i - 1] - mt[i])
            return t[i - 1] + f * (t[i] - t[i - 1])
    return None


def main(path, t_end):
    t, s, fr = series(path, t_end)
    fig = plt.figure(figsize=(13.5, 8.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.95], hspace=0.38, wspace=0.28)

    # --- A: CMT angle, the headline ------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    ax.axhspan(45, 90, color="#15706c", alpha=0.05)
    ax.axhline(45, color="#555", ls=":", lw=1)
    ax.plot(t, s["mt_in"], "-", color=C_IN, lw=2.2, label="inner flank")
    ax.plot(t, s["mt_out"], "-", color=C_OUT, lw=2.2, label="outer flank")
    for mt, col, nm in ((s["mt_in"], C_IN, "inner"), (s["mt_out"], C_OUT, "outer")):
        x = crossing(t, mt)
        if x is not None:
            ax.axvline(x, color=col, ls="--", lw=1, alpha=0.7)
            ax.annotate(f"{nm} switches\n{x:.1f} h", (x, 68 if nm == "outer" else 78),
                        xytext=(6, 0), textcoords="offset points",
                        fontsize=7.5, color=col, va="center")
    ax.set(xlabel="hours after illumination", ylim=(-3, 93), xlim=(0, t[-1]),
           ylabel="predicted CMT angle to organ axis (deg)",
           title="A  CMT reorientation (cf. Figs 3E–F)")
    ax.set_yticks([0, 45, 90])
    ax.set_yticklabels(["0\nlongitudinal", "45", "90\ncircumferential"])
    ax.legend(fontsize=8, loc="lower left")

    # --- B: signed stress anisotropy, the paper's own plot --------------
    ax = fig.add_subplot(gs[0, 1])
    ax.axhline(0, color="#555", lw=1)
    ax.plot(t, s["a_in"], "-", color=C_IN, lw=2.2, label="inner flank")
    ax.plot(t, s["a_out"], "-", color=C_OUT, lw=2.2, label="outer flank")
    ax.fill_between(t, 0, s["a_out"], where=s["a_out"] < 0, color=C_OUT, alpha=0.12)
    ax.set(xlabel="hours after illumination", xlim=(0, t[-1]),
           ylabel="signed stress anisotropy",
           title="B  Stress anisotropy (cf. Fig 4E)")
    ax.text(0.98, 0.95, "circumferential", transform=ax.transAxes, ha="right",
            va="top", fontsize=7.5, color="#555")
    ax.text(0.98, 0.04, "longitudinal", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=7.5, color="#555")
    ax.legend(fontsize=8, loc="center right")

    # --- C: opening kinetics against experiment -------------------------
    ax = fig.add_subplot(gs[0, 2])
    et, ea = ED.OPENING["t"], np.array(ED.OPENING["a"])
    ese = np.array(ED.OPENING["se"])
    ax.fill_between(et, ea - ese, ea + ese, color=C_EXP, alpha=0.15, lw=0)
    ax.plot(et, ea, "o-", color=C_EXP, ms=4, lw=1.4, label="experiment ± SE")
    ax.plot(t, s["angle"], "-", color=C_SIM, lw=2.2, label="simulation")
    def rmse_to(hi):
        return math.sqrt(np.mean([(np.interp(x, t, s["angle"]) - y) ** 2
                                  for x, y in zip(et, ea) if x <= hi]))
    ax.axvspan(0, 6, color=C_SIM, alpha=0.06, lw=0)
    ax.text(0.03, 0.12, f"RMSE {rmse_to(6):.1f}° to 6 h\n"
                        f"       {rmse_to(10):.1f}° to 10 h",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=8,
            color=C_SIM)
    ax.text(3.0, 172, "fitted range", fontsize=7.5, color=C_SIM, ha="center")
    ax.set(xlabel="hours after illumination", ylabel="hook angle (deg)",
           title="C  Opening kinetics", xlim=(0, max(10, t[-1])), ylim=(0, 180))
    ax.legend(fontsize=8, loc="lower left", bbox_to_anchor=(0.0, 0.22))

    # --- D-F: polar histograms before / at / after the switch -----------
    x_out = crossing(t, s["mt_out"]) or t[-1] * 0.5
    picks = [0.0, min(t[-1], x_out + 1.0), t[-1]]
    labels = ["dark (t = 0)", f"after switch ({picks[1]:.1f} h)",
              f"late ({picks[2]:.0f} h)"]
    for col, (tp, lab) in enumerate(zip(picks, labels)):
        ax = fig.add_subplot(gs[1, col], projection="polar")
        idx = int(np.argmin(np.abs(t - tp)))
        v, cells = fr[idx]
        for sec, colr, nm in ((INNER, C_IN, "inner"), (OUTER, C_OUT, "outer")):
            angs = []
            mid = AB.N_BASAL + AB.N_HOOK // 2
            for i in range(mid - 6, mid + 6):
                p, q = AB.ring_center(v, i), AB.ring_center(v, i + 1)
                axis = tuple(b - a for a, b in zip(p, q))
                for j in sec:
                    c = cells[i * NC + j]
                    nvec = c[AB.SDIR:AB.SDIR + 3]
                    if math.hypot(*nvec) < 1e-12:
                        continue
                    a = AB.ang(axis, nvec)
                    angs.append(min(a, 180 - a))
            h, edges = np.histogram(angs, bins=18, range=(0, 90))
            centres = np.deg2rad(0.5 * (edges[:-1] + edges[1:]))
            # offset the two flanks so neither hides the other
            off = np.deg2rad(-1.3 if nm == "inner" else 1.3)
            ax.bar(centres + off, h, width=np.deg2rad(2.4), color=colr,
                   alpha=0.85, edgecolor="none", label=nm)
        ax.set_thetamin(0); ax.set_thetamax(90)
        ax.set_xticks(np.deg2rad([0, 45, 90]))
        ax.set_xticklabels(["0°\nlong.", "45°", "90°\ncirc."], fontsize=7.5)
        ax.set_yticklabels([])
        ax.set_title(f"{'DEF'[col]}  {lab}", fontsize=9, pad=12)
        if col == 0:
            ax.legend(fontsize=7.5, loc="upper right", bbox_to_anchor=(1.25, 1.1))

    fig.suptitle("Sub-epidermal growth drives the CMT switch, outer flank first",
                 y=0.985, fontsize=12)
    fig.text(0.005, 0.004, "simulation: InnerTissue::GrowingCore at force balance "
             "(QuasiStatic, no drag lag).  experiment: Walia, Carter et al. 2024, "
             "OpeningData; CMT angles are qualitative in Figs 3E–F (no tabulated "
             "values published)", fontsize=6.5, color="#7b8493")
    fig.savefig(f"{OUT}/fig7_mt_switch.png", dpi=150, bbox_inches="tight")
    print(f"{OUT}/fig7_mt_switch.png")
    for nm, mt in (("inner", s["mt_in"]), ("outer", s["mt_out"])):
        x = crossing(t, mt)
        print(f"  {nm} switch: {'%.2f h' % x if x else 'not within the run'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "run_core.out",
         float(sys.argv[2]) if len(sys.argv) > 2 else 12.0)
