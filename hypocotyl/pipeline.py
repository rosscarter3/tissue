#!/usr/bin/env python3
"""Equilibrate the biological hook shell in the dark, then run light and dark."""
import os, subprocess, sys
SIM = '/Users/ross/projects/tissue/build/simulator'
BASE = 'hook.model'

import re

def _zero(m, tag):
    """Set the parameter value on the line carrying @tag to zero."""
    return re.sub(r'^[0-9.eE+-]+(\s*#\s*@' + tag + r')', r'0.0\1', m, flags=re.M)

def variant(out, nogrowth=False, dark=False):
    m = open(BASE).read()
    if nogrowth:
        m = _zero(m, 'GROWTH')
        m = _zero(m, 'GROWTH_CT')
    if dark:
        # No illumination: auxin is not depleted, and the apoplast neither
        # acidifies nor relaxes - it holds the measured dark state in which
        # the inner flank is more alkaline than the outer.
        m = _zero(m, 'AUXIN_DECAY')
        m = _zero(m, 'ACID_DRIVE')
        m = _zero(m, 'ACID_RELAX')
    open(out, 'w').write(m)

if __name__ == '__main__':
    solver = sys.argv[1] if len(sys.argv) > 1 else 'solver.rk5'
    open('solver_eq.rk5','w').write('RK5Adaptive\n0 2.0\n0 2\n0.02 1e-2\n')
    variant('equil.model', nogrowth=True, dark=True)
    subprocess.run([SIM, 'equil.model', 'hook.init', 'solver_eq.rk5',
                    '-init_output', 'hook_eq.init',
                    '-init_output_format', 'centerTriTissue'],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    variant('hook_dark.model', dark=True)
    jobs = [(BASE, 'run_light.out', solver, []),
            ('hook_dark.model', 'run_dark.out', solver, [])]
    if os.path.exists('solver_vtk.rk5'):
        jobs.append((BASE, os.devnull, 'solver_vtk.rk5',
                     ['-vtk_output', os.path.join(os.getcwd(), 'vtk')]))
    procs = [subprocess.Popen([SIM, m, 'hook_eq.init', s, '-centerTri_init'] + extra,
                              stdout=open(o, 'w'), stderr=subprocess.DEVNULL)
             for m, o, s, extra in jobs]
    for p in procs: p.wait()
    print("done")
