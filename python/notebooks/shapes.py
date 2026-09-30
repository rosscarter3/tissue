"""Two shape statistics, so the notebook can say something about its output.

Deliberately small and self-contained: the point of the notebook is the
simulator, not a shape library.
"""
import math


def _area_perimeter(loop):
    pts = [(float(x), float(y)) for x, y in loop]
    a = p = 0.0
    for k, (x1, y1) in enumerate(pts):
        x2, y2 = pts[(k + 1) % len(pts)]
        a += x1 * y2 - x2 * y1
        p += math.hypot(x2 - x1, y2 - y1)
    return abs(a) / 2.0, p


def circularity(loop):
    """4*pi*A/P^2: 1 for a disc, falling as the outline gets longer."""
    a, p = _area_perimeter(loop)
    return 4 * math.pi * a / p ** 2 if p else float("nan")


def _hull(pts):
    pts = sorted(set((round(x, 9), round(y, 9)) for x, y in pts))
    if len(pts) < 3:
        return pts

    def half(ps):
        out = []
        for q in ps:
            while len(out) >= 2:
                (x1, y1), (x2, y2) = out[-2], out[-1]
                if (x2 - x1) * (q[1] - y1) - (y2 - y1) * (q[0] - x1) > 0:
                    break
                out.pop()
            out.append(q)
        return out

    return half(pts)[:-1] + half(pts[::-1])[:-1]


def solidity(loop):
    """area / convex-hull area.

    Built from areas rather than a perimeter, so unlike circularity it barely
    moves when the outline is drawn at a different resolution.
    """
    a, _ = _area_perimeter(loop)
    h, _ = _area_perimeter(_hull([tuple(p) for p in loop]))
    return a / h if h else float("nan")
