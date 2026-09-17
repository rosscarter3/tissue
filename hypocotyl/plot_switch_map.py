#!/usr/bin/env python3
"""Where and when the CMT switch happens, cell by cell.

The analytic model of Walia, Carter et al. (2024) resolves the hook into a few
independent parts at chosen angles phi. A cell-resolved shell instead gives a
switch *time* for every cell, so the switch becomes a map over the organ rather
than three curves - and the map is a prediction that can be compared against
per-cell CMT imaging.

Writes figures/fig8_switch_map.png.
"""
import math, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import analyze as AB

NC = AB.NCIRC
OUT = "figures"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                     "figure.facecolor": "white", "legend.frameon": False})


def per_cell_angles(path, t_end):
    fr = AB.parse(path)
    _, n_print = AB.print_times(t_end)
    times = np.array([t_end * i / n_print for i in range(len(fr))])
    nrings = (len(fr[0][1]) - 2) // NC
    ang = np.full((len(fr), nrings, NC), np.nan)
    for k, (v, cells) in enumerate(fr):
        for i in range(nrings):
            p, q = AB.ring_center(v, i), AB.ring_center(v, i + 1)
            axis = tuple(b - a for a, b in zip(p, q))
            for j in range(NC):
                n = cells[i * NC + j][AB.SDIR:AB.SDIR + 3]
                if math.hypot(*n) < 1e-12:
                    continue
                a = AB.ang(axis, n)
                ang[k, i, j] = min(a, 180 - a)
    return times, ang


def switch_times(times, ang):
    nr, nc = ang.shape[1], ang.shape[2]
    out = np.full((nr, nc), np.nan)
    for i in range(nr):
        for j in range(nc):
            s = ang[:, i, j]
            for k in range(1, len(s)):
                if s[k - 1] > 45 >= s[k]:
                    f = (s[k - 1] - 45) / (s[k - 1] - s[k])
                    out[i, j] = times[k - 1] + f * (times[k] - times[k - 1])
                    break
    return out


def main(path="run_core.out", t_end=12.0):
    os.makedirs(OUT, exist_ok=True)
    times, ang = per_cell_angles(path, t_end)
    sw = switch_times(times, ang)
    hook = slice(AB.N_BASAL, AB.N_BASAL + AB.N_HOOK)

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.3),
                             gridspec_kw={"width_ratios": [1.25, 1, 1]})

    # --- map over the organ --------------------------------------------
    ax = axes[0]
    # circumferential index 0 = inner flank, NC/2 = outer flank
    phi = (np.arange(NC) / NC) * 360.0
    im = ax.imshow(sw[hook].T, origin="lower", aspect="auto", cmap="viridis",
                   extent=[AB.N_BASAL, AB.N_BASAL + AB.N_HOOK, 0, 360])
    ax.set(xlabel="ring index along the hook (base → apex)",
           ylabel="circumferential position (deg)",
           title="A  switch time, cell by cell")
    ax.set_yticks([0, 90, 180, 270, 360])
    ax.set_yticklabels(["0 inner", "90", "180 outer", "270", "360 inner"])
    cb = fig.colorbar(im, ax=ax); cb.set_label("hours to switch")
    nsw = np.isfinite(sw[hook]).sum()
    ax.text(0.02, 0.96, f"{nsw} of {sw[hook].size} hook cells switch",
            transform=ax.transAxes, va="top", fontsize=7.5, color="white")

    # --- switch time around the circumference --------------------------
    ax = axes[1]
    mid = AB.N_BASAL + AB.N_HOOK // 2
    band = sw[mid - 6:mid + 6]
    med = np.nanmedian(band, axis=0)
    lo = np.nanpercentile(band, 10, axis=0)
    hi = np.nanpercentile(band, 90, axis=0)
    ax.fill_between(phi, lo, hi, color="#46337e", alpha=0.18, lw=0)
    ax.plot(phi, med, "o-", color="#46337e", ms=4, lw=1.6)
    ax.axvline(0, color="#15706c", ls=":", lw=1.2)
    ax.axvline(180, color="#b04d22", ls=":", lw=1.2)
    ax.text(6, ax.get_ylim()[1], " inner", color="#15706c", va="top", fontsize=8)
    ax.text(186, ax.get_ylim()[1], " outer", color="#b04d22", va="top", fontsize=8)
    ax.set(xlabel="circumferential position (deg)", ylabel="hours to switch",
           title="B  switch time around the circumference", xlim=(0, 360))
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.grid(alpha=0.16)

    # --- the analytic prediction it is being compared with --------------
    ax = axes[2]
    lam = 2.88
    ph = np.linspace(0, 2 * np.pi, 361)
    hoop = (2 * lam + np.sin(ph - np.pi / 2)) / (lam + np.sin(ph - np.pi / 2))
    ax.plot(np.rad2deg(ph), hoop, "-", color="#22272e", lw=2,
            label=r"analytic hoop stress")
    ax.set(xlabel="circumferential position (deg)",
           ylabel="hoop stress / reference",
           title="C  what sets the ordering (analytic toroid)", xlim=(0, 360))
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.grid(alpha=0.16)
    ax.legend(fontsize=8)
    ax.axvline(180, color="#b04d22", ls=":", lw=1.2)
    ax.annotate("least hoop stress here,\nso overtaken first", xy=(180, 1.742),
                xytext=(180, 2.15), ha="center", fontsize=7.5, color="#b04d22",
                arrowprops=dict(arrowstyle="->", color="#b04d22", lw=1))

    fig.suptitle("A cell-resolved switch map: a prediction the analytic model "
                 "cannot make", y=1.0, fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(f"{OUT}/fig8_switch_map.png", dpi=150, bbox_inches="tight")
    print(f"{OUT}/fig8_switch_map.png")
    print(f"  hook cells switching: {nsw}/{sw[hook].size}")
    print(f"  earliest {np.nanmin(sw[hook]):.2f} h at phi="
          f"{phi[np.nanargmin(np.nanmin(sw[hook], axis=0))]:.0f} deg")
    print(f"  spread across the switching band: "
          f"{np.nanmax(sw[hook]) - np.nanmin(sw[hook]):.2f} h")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "run_core.out",
         float(sys.argv[2]) if len(sys.argv) > 2 else 12.0)
