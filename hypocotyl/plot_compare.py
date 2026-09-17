#!/usr/bin/env python3
"""Reproduce the experimental plots of Walia, Carter et al. (2024) Dev Cell
59:3245 and overlay the tissue-v2 simulation.

Experimental values come from the authors' published analysis package
(zenodo.org/records/13379829); simulation values are read from print-flag-0
output of build/simulator. Figures are written to ./figures/.
"""
import math, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import expdata as ED
import analyze as AB

OUT = "figures"
os.makedirs(OUT, exist_ok=True)

# one palette for the whole set: experiment dark, simulation coloured
C_EXP   = "#22272e"
C_SIM   = "#46337e"
C_INNER = "#15706c"
C_OUTER = "#b04d22"
C_DARK  = "#7b8493"
plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.grid": True,
    "grid.alpha": 0.16, "grid.linewidth": 0.6, "legend.frameon": False,
})

# ----------------------------------------------------------------- simulation
class Sim:
    """Per-timepoint observables from one print-flag-0 run."""
    def __init__(self, path, t_end=12.0):
        self.frames = AB.parse(path)
        n = len(self.frames)
        # time from the solver's print schedule, not the frame count, so
        # an unfinished run is not stretched over the whole interval
        _, n_print = AB.print_times(t_end)
        self.t = np.array([t_end * i / n_print for i in range(n)])
        self.complete = n >= n_print + 1
        self.angle, self.inner, self.outer, self.middle = [], [], [], []
        self.mt_in, self.mt_out, self.an_in, self.an_out = [], [], [], []
        self.aux_in, self.aux_out, self.acid_in, self.acid_out = [], [], [], []
        NC = AB.NCIRC
        inner_s, outer_s, mid_s = [0, NC - 1], [NC // 2 - 1, NC // 2], [NC // 4]
        a0 = AB.arc(self.frames[0][0], 0)
        o0 = AB.arc(self.frames[0][0], NC // 2)
        m0 = AB.arc(self.frames[0][0], NC // 4)
        for v, c in self.frames:
            self.angle.append(AB.hook_angle(v))
            self.inner.append(AB.arc(v, 0) / a0)
            self.outer.append(AB.arc(v, NC // 2) / o0)
            self.middle.append(AB.arc(v, NC // 4) / m0)
            mi, ai, xi, ci = AB.flank_mt(v, c, inner_s)
            mo, ao, xo, co = AB.flank_mt(v, c, outer_s)
            self.mt_in.append(mi);  self.mt_out.append(mo)
            self.an_in.append(ai);  self.an_out.append(ao)
            self.aux_in.append(xi); self.aux_out.append(xo)
            self.acid_in.append(ci); self.acid_out.append(co)
        for k in ("angle inner outer middle mt_in mt_out an_in an_out "
                  "aux_in aux_out acid_in acid_out").split():
            setattr(self, k, np.array(getattr(self, k)))

    def flank_angles(self, frame_idx, sectors):
        """Per-cell CMT angles at one timepoint, for polar histograms."""
        v, cells = self.frames[frame_idx]
        out, mid = [], AB.N_BASAL + AB.N_HOOK // 2
        for i in range(mid - 6, mid + 6):
            p, q = AB.ring_center(v, i), AB.ring_center(v, i + 1)
            ax = tuple(b - a for a, b in zip(p, q))
            for j in sectors:
                c = cells[i * AB.NCIRC + j]
                n = c[AB.SDIR:AB.SDIR + 3]
                if math.hypot(*n) < 1e-12:
                    continue
                a = AB.ang(ax, n)
                out.append(min(a, 180 - a))
        return np.array(out)

def load(path, t_end=12.0):
    return Sim(path, t_end) if os.path.exists(path) else None

def stamp(fig, text):
    fig.text(0.005, 0.004, text, fontsize=6.5, color="#7b8493", ha="left")

# ------------------------------------------------------------------ figure 1
def fig_opening(light, dark):
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0))
    ax = axes[0]
    t, a, se = ED.OPENING["t"], np.array(ED.OPENING["a"]), np.array(ED.OPENING["se"])
    ax.fill_between(t, a - se, a + se, color=C_EXP, alpha=0.15, lw=0)
    ax.plot(t, a, "o-", color=C_EXP, ms=4, lw=1.4, label="experiment (mean ± SE)")
    if light is not None:
        ax.plot(light.t, light.angle, "-", color=C_SIM, lw=2, label="simulation, light")
    if dark is not None:
        ax.plot(dark.t, dark.angle, "--", color=C_DARK, lw=1.8, label="simulation, dark")
    ax.set(xlabel="hours after illumination", ylabel="hook angle (deg)",
           title="Hook opening kinetics", xlim=(0, 10), ylim=(0, 180))
    ax.legend(fontsize=8, loc="upper right")
    if light is not None:
        rmse = math.sqrt(np.mean([(np.interp(ti, light.t, light.angle) - ai) ** 2
                                  for ti, ai in zip(t, a)]))
        ax.text(0.03, 0.05, f"RMSE {rmse:.1f}°", transform=ax.transAxes,
                ha="left", va="bottom", fontsize=8, color=C_SIM)

    ax = axes[1]
    lt = ED.LIGHT["t"]
    for key, col, lab in (("normal_light", C_EXP, "normal light"),
                          ("low_light", "#8a7fb8", "low light"),
                          ("dark", C_DARK, "dark")):
        ax.plot(*ED.series(ED.LIGHT, key), "-", color=col, lw=1.5,
                label=f"exp, {lab}")
    if light is not None:
        ax.plot(light.t, light.angle, "-", color=C_SIM, lw=2, label="sim, light")
    if dark is not None:
        ax.plot(dark.t, dark.angle, "--", color=C_SIM, lw=1.6, alpha=0.7, label="sim, dark")
    ax.set(xlabel="hours after illumination", ylabel="hook angle (deg)",
           title="Light dose dependence", xlim=(0, 28), ylim=(0, 180))
    ax.legend(fontsize=7.5, ncol=2)
    fig.suptitle("Hook opening: experiment vs simulation", y=0.99, fontsize=11)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    stamp(fig, "experiment: Walia, Carter et al. 2024 Dev Cell 59:3245 (OpeningData, LightData)")
    fig.savefig(f"{OUT}/fig1_opening.png", dpi=150)
    print("fig1_opening.png")

# ------------------------------------------------------------------ figure 2
def fig_lengths(light):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    ax = axes[0]
    th = ED.FOLD["t_h"]
    for key, std, col, lab in (("inner_average", "inner_std", C_INNER, "inner"),
                               ("middle_average", "middle_std", "#6a7a86", "middle"),
                               ("outer_average", "outer_std", C_OUTER, "outer")):
        m, s = np.array(ED.FOLD[key]), np.array(ED.FOLD[std])
        ax.errorbar(th, m, yerr=s, fmt="o-", color=col, ms=4, lw=1.4,
                    capsize=2.5, label=f"exp {lab}")
    if light is not None:
        ax.plot(light.t, light.inner, "-", color=C_INNER, lw=2, alpha=0.55)
        ax.plot(light.t, light.middle, "-", color="#6a7a86", lw=2, alpha=0.55)
        ax.plot(light.t, light.outer, "-", color=C_OUTER, lw=2, alpha=0.55)
        ax.plot([], [], "-", color="k", alpha=0.55, lw=2, label="simulation")
    ax.set(xlabel="hours after illumination", ylabel="cell length fold change",
           title="Cell length by flank (Fig 1)", xlim=(0, 10))
    ax.legend(fontsize=7.5, ncol=2)
    ax.text(0.5, 0.03, "experiment is single-cell length; the simulation curve\n"
                       "is a tissue-level arc measure (middle panel)",
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=6.8, color="#7b8493")

    ax = axes[1]
    L = ED.lengths_csv()
    ax.plot(L["t"], L["inner"], "o-", color=C_INNER, ms=3.5, lw=1.3, label="exp inner")
    ax.plot(L["t"], L["outer_short"], "o-", color=C_OUTER, ms=3.5, lw=1.3, label="exp outer")
    if light is not None:
        ax.plot(light.t, light.inner, "-", color=C_INNER, lw=2, alpha=0.55, label="sim inner")
        ax.plot(light.t, light.outer, "-", color=C_OUTER, lw=2, alpha=0.55, label="sim outer")
    ax.set(xlabel="hours after illumination", ylabel="tissue length fold change",
           title="Tissue length time course", xlim=(0, 10))
    ax.legend(fontsize=7.5)

    ax = axes[2]
    rt = [t for t, r in zip(L["t"], L["inner_rate"]) if r is not None]
    ri = [r for r in L["inner_rate"] if r is not None]
    ro = [r for r in L["outer_rate"] if r is not None]
    ax.plot(rt, ri, "o-", color=C_INNER, ms=3.5, lw=1.2, label="exp inner")
    ax.plot(rt[:len(ro)], ro, "o-", color=C_OUTER, ms=3.5, lw=1.2, label="exp outer")
    if light is not None:
        # convert fold change to the same um/h units as the measurement
        li0, lo0 = 236.251, 572.968
        ax.plot(light.t[1:], np.diff(light.inner * li0) / np.diff(light.t),
                "-", color=C_INNER, lw=2, alpha=0.55, label="sim inner")
        ax.plot(light.t[1:], np.diff(light.outer * lo0) / np.diff(light.t),
                "-", color=C_OUTER, lw=2, alpha=0.55, label="sim outer")
    ax.axhline(0, color="#aab", lw=0.8)
    ax.set(xlabel="hours after illumination", ylabel="growth rate (um/h)",
           title="Growth rate", xlim=(0, 10))
    ax.legend(fontsize=7.5)
    fig.suptitle("Differential growth: experiment vs simulation", y=0.99, fontsize=11)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    stamp(fig, "experiment: LengthFoldChangeData and length_measurements.csv (Walia et al. 2024)")
    fig.savefig(f"{OUT}/fig2_lengths.png", dpi=150)
    print("fig2_lengths.png")

# ------------------------------------------------------------------ figure 3
def fig_chemistry(light, dark):
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0))
    ax = axes[0]
    if light is not None:
        ax.plot(light.t, light.aux_in, "-", color=C_INNER, lw=2, label="sim inner")
        ax.plot(light.t, light.aux_out, "-", color=C_OUTER, lw=2, label="sim outer")
    if dark is not None:
        ax.plot(dark.t, dark.aux_in, "--", color=C_INNER, lw=1.4, alpha=0.7, label="dark inner")
        ax.plot(dark.t, dark.aux_out, "--", color=C_OUTER, lw=1.4, alpha=0.7, label="dark outer")
    ax.set(xlabel="hours after illumination", ylabel="auxin (relative)",
           title="Auxin: DR5 maximum inner, depleted by light", xlim=(0, 10), ylim=(0, 1.05))
    ax.annotate("measured: inner DR5 max in dark,\ndepleted after illumination",
                xy=(0.97, 0.72), xycoords="axes fraction", ha="right",
                fontsize=7.5, color=C_EXP)
    ax.legend(fontsize=7.5, ncol=2)

    ax = axes[1]
    if light is not None:
        ax.plot(light.t, light.acid_in, "-", color=C_INNER, lw=2, label="sim inner")
        ax.plot(light.t, light.acid_out, "-", color=C_OUTER, lw=2, label="sim outer")
    if dark is not None:
        ax.plot(dark.t, dark.acid_in, "--", color=C_INNER, lw=1.4, alpha=0.7, label="dark inner")
        ax.plot(dark.t, dark.acid_out, "--", color=C_OUTER, lw=1.4, alpha=0.7, label="dark outer")
    ax.axvline(0.5, color=C_EXP, lw=1, ls=":")
    ax.annotate("measured: both flanks\nacidify by 30 min", xy=(0.55, 0.30),
                xytext=(2.4, 0.30), fontsize=7.5, color=C_EXP,
                arrowprops=dict(arrowstyle="->", color=C_EXP, lw=0.9))
    ax.set(xlabel="hours after illumination", ylabel="apoplastic acidification (0 = dark state)",
           title="Apoplastic pH: inner more alkaline in dark", xlim=(0, 10), ylim=(0, 1.05))
    ax.legend(fontsize=7.5, ncol=2)
    fig.suptitle("Chemical gates (Fig 2): both released by light", y=0.99, fontsize=11)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    stamp(fig, "the paper reports pH as a relative biosensor ratio, so the measured landmarks are annotated rather than overlaid")
    fig.savefig(f"{OUT}/fig3_chemistry.png", dpi=150)
    print("fig3_chemistry.png")

# ------------------------------------------------------------------ figure 4
def fig_anisotropy(light, dark):
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0))
    ax = axes[0]
    if light is not None:
        ax.plot(light.t, light.mt_in, "-", color=C_INNER, lw=2, label="sim inner")
        ax.plot(light.t, light.mt_out, "-", color=C_OUTER, lw=2, label="sim outer")
    if dark is not None:
        ax.plot(dark.t, dark.mt_in, "--", color=C_INNER, lw=1.4, alpha=0.7, label="dark inner")
        ax.plot(dark.t, dark.mt_out, "--", color=C_OUTER, lw=1.4, alpha=0.7, label="dark outer")
    ax.axhline(45, color=C_EXP, lw=1, ls=":")
    ax.text(9.7, 47, "transverse / longitudinal boundary", ha="right",
            fontsize=7, color=C_EXP)
    ax.set(xlabel="hours after illumination",
           ylabel="CMT angle to organ axis (deg)",
           title="Predicted CMT orientation (Fig 3E-F)", xlim=(0, 10), ylim=(0, 95))
    ax.set_yticks([0, 45, 90])
    ax.set_yticklabels(["0 longitudinal", "45", "90 circumferential"])
    if light is not None:
        below = np.where(light.mt_out < 45)[0]
        if len(below):
            ax.annotate(f"outer switches at {light.t[below[0]]:.1f} h",
                        xy=(light.t[below[0]], 45), xytext=(light.t[below[0]] + 1.1, 63),
                        fontsize=7.5, color=C_OUTER,
                        arrowprops=dict(arrowstyle="->", color=C_OUTER, lw=0.9))
    ax.legend(fontsize=7.5, ncol=2)

    ax = axes[1]
    if light is not None:
        ax.plot(light.t, light.an_in, "-", color=C_INNER, lw=2, label="sim inner")
        ax.plot(light.t, light.an_out, "-", color=C_OUTER, lw=2, label="sim outer")
    if dark is not None:
        ax.plot(dark.t, dark.an_in, "--", color=C_INNER, lw=1.4, alpha=0.7, label="dark inner")
        ax.plot(dark.t, dark.an_out, "--", color=C_OUTER, lw=1.4, alpha=0.7, label="dark outer")
    ax.set(xlabel="hours after illumination", ylabel="stress anisotropy  a = 1 - s2/s1",
           title="Stress anisotropy (Fig 3D, 4E)", xlim=(0, 10), ylim=(0, 1))
    ax.legend(fontsize=7.5, ncol=2)
    fig.suptitle("Microtubule orientation and stress anisotropy", y=0.99, fontsize=11)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    stamp(fig, "paper: CMTs circumferential in dark, switching to longitudinal after light, outer flank before inner")
    fig.savefig(f"{OUT}/fig4_anisotropy.png", dpi=150)
    print("fig4_anisotropy.png")

# ------------------------------------------------------------------ figure 5
def fig_polar(light):
    if light is None:
        return
    n = len(light.frames)
    picks = [(0, "0 h"), (n // 6, f"{light.t[n//6]:.0f} h"),
             (n // 3, f"{light.t[n//3]:.0f} h"), (n - 1, f"{light.t[-1]:.0f} h")]
    fig, axes = plt.subplots(2, len(picks), figsize=(3.0 * len(picks), 6.2),
                             subplot_kw={"projection": "polar"})
    NC = AB.NCIRC
    for row, (sectors, lab, col) in enumerate(
            (([0, NC - 1], "inner", C_INNER), ([NC // 2 - 1, NC // 2], "outer", C_OUTER))):
        for col_i, (fi, tlab) in enumerate(picks):
            ax = axes[row, col_i]
            angs = light.flank_angles(fi, sectors)
            bins = np.linspace(0, 90, 10)
            h, edges = np.histogram(angs, bins=bins)
            widths = np.diff(np.radians(edges))
            ax.bar(np.radians(edges[:-1]), h, width=widths, align="edge",
                   color=col, alpha=0.85, edgecolor="white", linewidth=0.5)
            ax.set_thetamin(0); ax.set_thetamax(90)
            ax.set_xticks(np.radians([0, 45, 90]))
            ax.set_xticklabels(["0 long.", "45", "90 circ."], fontsize=7)
            ax.set_yticklabels([])
            ax.set_title(f"{lab}, {tlab}", fontsize=9, color=col, pad=12)
    fig.suptitle("Predicted CMT angle distributions (cf. Fig 3E polar histograms)",
                 y=0.98, fontsize=11)
    fig.tight_layout(rect=(0, 0.02, 1, 0.95))
    stamp(fig, "simulation only: per-cell maximal principal stress direction, mid-hook cells")
    fig.savefig(f"{OUT}/fig5_cmt_polar.png", dpi=150)
    print("fig5_cmt_polar.png")

# ------------------------------------------------------------------ figure 6
def fig_perturbations(light):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    t = ED.ISOX["t"]
    ax = axes[0]
    doses = (("mock", C_EXP, "mock"), ("isox_100nm", "#c77", "isoxaben 100 nM"),
             ("isox_600nm", "#944", "isoxaben 600 nM"))
    for key, col, lab in doses:
        ax.plot(*ED.series(ED.ISOX, key), "o-", color=col, ms=3.5, lw=1.3,
                label=f"exp {lab}")
    for name, col, lab in (("pt_isox100.out", "#c77", "sim isox 100 nM"),
                           ("pt_isox600.out", "#944", "sim isox 600 nM")):
        s = load(name)
        if s: ax.plot(s.t, s.angle, "-", color=col, lw=2, alpha=0.55, label=lab)
    if light is not None:
        ax.plot(light.t, light.angle, "-", color=C_SIM, lw=2, alpha=0.6, label="sim mock")
    ax.set(xlabel="hours after illumination", ylabel="hook angle (deg)",
           title="Isoxaben: less cellulose (lower Y_fiber)", xlim=(0, 10), ylim=(0, 180))
    ax.legend(fontsize=7)
    ax.text(0.5, 0.06, "model: no effect — the fibre carries too little\nload for wall anisotropy to matter (see README)",
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=6.8, color="#a33", style="italic")

    ax = axes[1]
    for key, col, lab in (("mock", C_EXP, "mock"), ("oryzalin", "#7a5", "oryzalin")):
        ax.plot(*ED.series(ED.ORYZALIN, key), "o-", color=col, ms=3.5, lw=1.3,
                label=f"exp {lab}")
    s = load("pt_oryzalin.out")
    if s: ax.plot(s.t, s.angle, "-", color="#7a5", lw=2, alpha=0.6, label="sim oryzalin")
    if light is not None:
        ax.plot(light.t, light.angle, "-", color=C_SIM, lw=2, alpha=0.6, label="sim mock")
    ax.set(xlabel="hours after illumination", ylabel="hook angle (deg)",
           title="Oryzalin: no CMT guidance (isotropic fibre)", xlim=(0, 10), ylim=(0, 180))
    ax.legend(fontsize=7)
    ax.text(0.5, 0.06, "model: no effect — the fibre carries too little\nload for wall anisotropy to matter (see README)",
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=6.8, color="#a33", style="italic")

    ax = axes[2]
    Y = ED.YUC6
    if Y:
        ax.plot(Y["t"], Y["col0_a"], "o-", color=C_EXP, ms=3.5, lw=1.3, label="exp Col-0")
        ax.plot(Y["t"], Y["yuc6_a"], "o-", color="#b58", ms=3.5, lw=1.3, label="exp YUC6-OX")
    s = load("pt_yuc6.out")
    if s: ax.plot(s.t, s.angle, "-", color="#b58", lw=2, alpha=0.6, label="sim YUC6-OX")
    if light is not None:
        ax.plot(light.t, light.angle, "-", color=C_SIM, lw=2, alpha=0.6, label="sim WT")
    ax.set(xlabel="hours after illumination", ylabel="hook angle (deg)",
           title="YUC6-OX: auxin gate never released", xlim=(0, 10), ylim=(0, 180))
    ax.legend(fontsize=7)
    fig.suptitle("Perturbations: each treatment changes one model parameter "
             "(auxin gates reproduce; wall-fibre treatments do not)", y=0.99, fontsize=11)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    stamp(fig, "experiment: IsoxabenData, OryzalinData, YUC6Data (Walia et al. 2024)")
    fig.savefig(f"{OUT}/fig6_perturbations.png", dpi=150)
    print("fig6_perturbations.png")

if __name__ == "__main__":
    T = float(sys.argv[1]) if len(sys.argv) > 1 else 12.0
    light = load("run_light.out", T)
    dark = load("run_dark.out", T)
    if light is None:
        sys.exit("run_light.out not found - run pipeline.py first")
    for name in ("run_light.out", "run_dark.out") + tuple(
            f"pt_{k}.out" for k in ("oryzalin", "isox100", "isox600",
                                    "yuc6", "lowlight")):
        s = load(name, T)
        if s is not None and not s.complete:
            print(f"  note: {name} is still running "
                  f"(reaches t={s.t[-1]:.2f} of {T:g} h)")
    fig_opening(light, dark)
    fig_lengths(light)
    fig_chemistry(light, dark)
    fig_anisotropy(light, dark)
    fig_polar(light)
    fig_perturbations(light)
    print(f"\nfigures written to ./{OUT}/")
