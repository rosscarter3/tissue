#!/usr/bin/env python3
"""Equilibrate the 3D hook shell (dark, growth frozen), then light/dark runs."""
import subprocess, sys
SIM = '/Users/ross/projects/tissue/v2/build/simulator'

def variant(src, out, nogrowth=False, dark=False):
    m = open(src).read()
    if nogrowth:
        m = m.replace('20.0    # k_growth (1/h per unit strain): Lockhart wall yielding',
                      '0.0    # k_growth: frozen for equilibration')
        m = m.replace('20.0    # k_growth: internal edges follow the same yielding law',
                      '0.0    # k_growth: frozen for equilibration')
    if dark:
        m = m.replace('0.9     # k_d: light-triggered auxin depletion (t=0 = illumination)',
                      '0.0     # k_d: dark')
    open(out, 'w').write(m)

if __name__ == '__main__':
    solver = sys.argv[1] if len(sys.argv) > 1 else 'solver3d.rk5'
    variant('hook3d.model', 'equil3d.model', nogrowth=True, dark=True)
    open('solver3d_eq.rk5','w').write('RK5Adaptive\n0 2.5\n0 2\n0.02 1e-4\n')
    subprocess.run([SIM, 'equil3d.model', 'hook3d.init', 'solver3d_eq.rk5',
                    '-init_output', 'hook3d_eq.init',
                    '-init_output_format', 'centerTriTissue'],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    variant('hook3d.model', 'hook3d_dark.model', dark=True)
    procs = []
    for model, out in (('hook3d.model','run3d_light.out'),
                       ('hook3d_dark.model','run3d_dark.out')):
        procs.append(subprocess.Popen([SIM, model, 'hook3d_eq.init', solver,
                                       '-centerTri_init'],
                     stdout=open(out,'w'), stderr=subprocess.DEVNULL))
    for p in procs:
        p.wait()
    print("done")
