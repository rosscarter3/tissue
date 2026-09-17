#!/usr/bin/env python3
"""Compare the biological hook run against the paper's measurements:
hook angle kinetics, inner/outer cell-length fold change, and the predicted
CMT (maximal principal stress) orientation on each flank."""
import math, sys

NCIRC = 16
N_BASAL, N_HOOK, N_APICAL = 3, 26, 3
NRINGS = N_BASAL + N_HOOK + N_APICAL + 1
SDIR, ANISO, AUXIN, ACID = 0, 3, 5, 6

# Walia et al. 2024, Data S1 (hook_model/data.py OpeningData) and
# length_measurements.csv
EXP_T = list(range(11))
EXP_ANGLE = [158.66,155.59,133.44,102.88,81.49,63.21,53.33,46.81,44.31,40.38,34.50]
EXP_INNER = {0:1.0,2:1.363,4:1.918,6:2.044,8:2.149}      # relative_inner
EXP_OUTER = {0:1.0,2:1.067,4:1.040,6:1.061,8:1.055}      # relative_outer_short

def print_times(t_end=12.0, solver="solver.rk5"):
    """Times of the print frames, from the solver's print schedule.

    Frame k is at t_end*k/n_print whatever the file length, so a run that is
    still in progress (or was cut short) is timed correctly instead of being
    stretched over the whole interval."""
    n_print = 48
    try:
        rows = [l.split() for l in open(solver) if l.strip()]
        t_end, n_print = float(rows[1][1]), int(rows[2][1])
    except (OSError, IndexError, ValueError):
        pass
    return t_end, n_print

def vid(i,j): return i*NCIRC + (j % NCIRC)

def parse(path):
    tok = open(path).read().split(); pos = 1; frames = []
    try:
        while pos < len(tok):
            nv, dim = int(tok[pos]), int(tok[pos+1]); pos += 2
            co = [float(x) for x in tok[pos:pos+nv*dim]]
            if len(co) < nv*dim: break
            pos += nv*dim
            verts = [tuple(co[dim*i:dim*i+3]) for i in range(nv)]
            nc, ncol = int(tok[pos]), int(tok[pos+1]); pos += 2
            cells = []
            for _ in range(nc):
                nvert = int(tok[pos])
                cells.append([float(x) for x in tok[pos+1+nvert:pos+1+nvert+ncol]])
                pos += 1 + nvert + ncol
            nw, wcol = int(tok[pos]), int(tok[pos+1]); pos += 2 + nw*(wcol+2)
            frames.append((verts, cells))
    except (ValueError, IndexError): pass
    return frames

def ring_center(v,i):
    p=[v[vid(i,j)] for j in range(NCIRC)]
    return tuple(sum(q[d] for q in p)/NCIRC for d in range(3))

def ang(u,w):
    d=sum(a*b for a,b in zip(u,w)); nu=math.hypot(*u); nw=math.hypot(*w)
    return math.degrees(math.acos(max(-1,min(1,d/(nu*nw)))))

def hook_angle(v):
    b=tuple(q-p for p,q in zip(ring_center(v,2),ring_center(v,6)))
    a=tuple(q-p for p,q in zip(ring_center(v,NRINGS-3),ring_center(v,NRINGS-1)))
    return ang(b,a)

def arc(v,j):
    return sum(math.dist(v[vid(i,j)],v[vid(i+1,j)]) for i in range(N_BASAL,N_BASAL+N_HOOK))

def flank_mt(v,cells,sectors):
    """Mean CMT angle to the organ axis: 90 deg = circumferential (transverse),
    0 deg = longitudinal. Also returns mean anisotropy, auxin, acid."""
    angs=[];ani=[];aux=[];aci=[]
    mid=N_BASAL+N_HOOK//2
    for i in range(mid-6,mid+6):
        p,q=ring_center(v,i),ring_center(v,i+1)
        ax=tuple(b-a for a,b in zip(p,q))
        for j in sectors:
            c=cells[i*NCIRC+j]
            n=c[SDIR:SDIR+3]
            if math.hypot(*n)<1e-12: continue
            a=ang(ax,n)
            angs.append(min(a,180-a))       # axis, not vector
            ani.append(c[ANISO]); aux.append(c[AUXIN]); aci.append(c[ACID])
    m=lambda x: sum(x)/len(x) if x else float('nan')
    return m(angs), m(ani), m(aux), m(aci)

if __name__=='__main__':
    path=sys.argv[1] if len(sys.argv)>1 else 'run_light.out'
    T=float(sys.argv[2]) if len(sys.argv)>2 else 12.0
    fr=parse(path); n=len(fr)
    inner=[0,NCIRC-1]; outer=[NCIRC//2-1,NCIRC//2]
    Tsolver,nprint=print_times(T)
    if len(sys.argv)<=2: T=Tsolver
    times=[T*k/nprint for k in range(n)]
    if n<nprint+1:
        print(f"# WARNING: {path} has {n}/{nprint+1} frames "
              f"(run reached t={times[-1]:.2f} of {T:g} h)\n")
    a0=arc(fr[0][0],0); o0=arc(fr[0][0],NCIRC//2)
    print(f"{'t':>5} {'angle':>7} {'exp':>6} | {'inner':>6} {'exp':>5} {'outer':>6} {'exp':>5} |"
          f" {'MTin':>5} {'MTout':>6} | {'auxin':>5} {'acid':>5}")
    for t in EXP_T:
        k=min(range(n),key=lambda x:abs(times[x]-t))
        v,c=fr[k]
        mi,ai,xi,ci=flank_mt(v,c,inner); mo,ao,xo,co=flank_mt(v,c,outer)
        ei=f"{EXP_INNER[t]:.2f}" if t in EXP_INNER else "  -"
        eo=f"{EXP_OUTER[t]:.2f}" if t in EXP_OUTER else "  -"
        print(f"{t:5d} {hook_angle(v):7.1f} {EXP_ANGLE[t]:6.1f} | {arc(v,0)/a0:6.2f} {ei:>5}"
              f" {arc(v,NCIRC//2)/o0:6.2f} {eo:>5} | {mi:5.0f} {mo:6.0f} | {xi:5.2f} {ci:5.2f}")
    rmse=math.sqrt(sum((hook_angle(fr[min(range(n),key=lambda x:abs(times[x]-t))][0])-EXP_ANGLE[t])**2
                       for t in EXP_T)/len(EXP_T))
    print(f"\nhook-angle RMSE vs experiment: {rmse:.1f} deg")
    print("MT angle: 90 = circumferential (transverse), 0 = longitudinal")
