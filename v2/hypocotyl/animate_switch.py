#!/usr/bin/env python3
"""Animate the hook opening with the predicted CMT switch visible.

Every epidermal cell is coloured by the angle between the maximal principal
stress direction - the axis the model predicts cortical microtubules follow -
and the local organ axis. Circumferential (transverse) cells are teal,
longitudinal cells are orange, so the outer flank visibly turns over partway
through the run while the inner flank does not.

Writes figures/hook_switch.mp4 (and a .gif fallback if ffmpeg is missing).

Usage: python3 animate_switch.py [run.out] [t_end_hours]
"""
import math, os, subprocess, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import analyze as AB

NC = AB.NCIRC
OUT = "figures"


def parse_with_faces(path):
    """Like analyze.parse but keeping each cell's vertex ring, for drawing."""
    tok = open(path).read().split()
    pos, frames = 1, []
    try:
        while pos < len(tok):
            nv, dim = int(tok[pos]), int(tok[pos + 1]); pos += 2
            co = [float(x) for x in tok[pos:pos + nv * dim]]
            if len(co) < nv * dim:
                break
            pos += nv * dim
            verts = np.array([co[dim * i:dim * i + 3] for i in range(nv)])
            nc, ncol = int(tok[pos]), int(tok[pos + 1]); pos += 2
            faces, cvars = [], []
            for _ in range(nc):
                n = int(tok[pos])
                faces.append([int(x) for x in tok[pos + 1:pos + 1 + n]])
                cvars.append([float(x) for x in tok[pos + 1 + n:pos + 1 + n + ncol]])
                pos += 1 + n + ncol
            nw, wcol = int(tok[pos]), int(tok[pos + 1]); pos += 2 + nw * (wcol + 2)
            frames.append((verts, faces, cvars))
    except (ValueError, IndexError):
        pass
    return frames


def cmt_angles(verts, faces, cvars):
    """Angle of each cell's stress direction to the local organ axis, 0-90."""
    out = np.full(len(faces), np.nan)   # caps stay NaN and are not drawn
    nrings = (len(faces) - 2) // NC
    for i in range(nrings):
        p = AB.ring_center(verts, i)
        q = AB.ring_center(verts, i + 1)
        axis = tuple(b - a for a, b in zip(p, q))
        for j in range(NC):
            c = i * NC + j
            n = cvars[c][AB.SDIR:AB.SDIR + 3]
            if math.hypot(*n) < 1e-12:
                continue
            a = AB.ang(axis, n)
            out[c] = min(a, 180 - a)
    return out


def main(path="run_core.out", t_end=12.0):
    os.makedirs(OUT, exist_ok=True)
    frames = parse_with_faces(path)
    n = len(frames)
    _, n_print = AB.print_times(t_end)
    times = [t_end * i / n_print for i in range(n)]
    print(f"{n} frames from {path}")

    cmap = LinearSegmentedColormap.from_list(
        "cmt", ["#c2521c", "#e0a050", "#dddddd", "#4aa39e", "#10605c"])

    # Global bounds across the whole run so the organ does not jump between
    # frames, with the box aspect matched to the data so it fills the figure.
    allv = np.vstack([f[0] for f in frames])
    lo = allv.min(axis=0) - 12.0
    hi = allv.max(axis=0) + 12.0
    ext = hi - lo

    fig = plt.figure(figsize=(7.6, 7.2), facecolor="#111114")
    ax = fig.add_subplot(111, projection="3d", facecolor="#111114")

    def draw(k):
        ax.clear()
        ax.set_facecolor("#111114")
        verts, faces, cvars = frames[k]
        ang = cmt_angles(verts, faces, cvars)
        polys, cols = [], []
        for c, f in enumerate(faces):
            if np.isnan(ang[c]):
                continue
            polys.append(verts[f])
            cols.append(cmap(ang[c] / 90.0))
        pc = Poly3DCollection(polys, facecolors=cols, edgecolors="#33333a",
                              linewidths=0.25)
        ax.add_collection3d(pc)
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
        ax.set_box_aspect(tuple(ext))
        ax.set_axis_off()
        ax.view_init(elev=14, azim=-72)
        ax.set_title(f"{times[k]:4.1f} h after illumination",
                     color="white", fontsize=15, pad=-6)
        inner = np.nanmean([ang[i * NC + j] for i in range(AB.N_BASAL + 6,
                            AB.N_BASAL + AB.N_HOOK - 6) for j in (0, NC - 1)])
        outer = np.nanmean([ang[i * NC + j] for i in range(AB.N_BASAL + 6,
                            AB.N_BASAL + AB.N_HOOK - 6) for j in (NC // 2 - 1, NC // 2)])
        ax.text2D(0.02, 0.10, f"inner flank  {inner:4.0f}°", transform=ax.transAxes,
                  color="#4aa39e", fontsize=11)
        ax.text2D(0.02, 0.05, f"outer flank  {outer:4.0f}°", transform=ax.transAxes,
                  color="#c2521c" if outer < 45 else "#4aa39e", fontsize=11)
        ax.text2D(0.62, 0.05, "teal = circumferential\norange = longitudinal",
                  transform=ax.transAxes, color="#8a8a92", fontsize=8.5)

    tmp = os.path.join(OUT, "_frames")
    os.makedirs(tmp, exist_ok=True)
    for k in range(n):
        draw(k)
        fig.savefig(os.path.join(tmp, f"f{k:04d}.png"), dpi=110,
                    facecolor="#111114")
    plt.close(fig)
    print(f"rendered {n} frames")

    mp4 = os.path.join(OUT, "hook_switch.mp4")
    try:
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "6",
                        "-i", os.path.join(tmp, "f%04d.png"),
                        "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                        "-c:v", "libx264", "-pix_fmt", "yuv420p", mp4], check=True)
        print(f"wrote {mp4}")
        for f in os.listdir(tmp):
            os.remove(os.path.join(tmp, f))
        os.rmdir(tmp)
    except (OSError, subprocess.CalledProcessError) as e:
        print(f"ffmpeg unusable ({type(e).__name__}); writing a GIF instead")
        from PIL import Image
        imgs = [Image.open(os.path.join(tmp, f"f{k:04d}.png")).convert("P",
                palette=Image.ADAPTIVE) for k in range(n)]
        gif = os.path.join(OUT, "hook_switch.gif")
        imgs[0].save(gif, save_all=True, append_images=imgs[1:], duration=170,
                     loop=0, optimize=True)
        print(f"wrote {gif}")
        for f in os.listdir(tmp):
            os.remove(os.path.join(tmp, f))
        os.rmdir(tmp)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "run_core.out",
         float(sys.argv[2]) if len(sys.argv) > 2 else 12.0)
