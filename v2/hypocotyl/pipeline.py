#!/usr/bin/env python3
"""Equilibrate the hook under turgor (dark, growth frozen), then run
light/dark from the equilibrated state."""
import subprocess, sys

SIM = '/Users/ross/projects/tissue/v2/build/simulator'

def variant(src, out, **kw):
    m = open(src).read()
    if kw.get('nogrowth'):
        m = m.replace('0.25    # k_growth (1/h): active elongation once auxin clears (acid growth;',
                      '0.0    # k_growth (1/h): active elongation once auxin clears (acid growth;')
        m = m.replace('0.15    # k_growth (1/h): active relative elongation rate',
                      '0.0    # k_growth (1/h): active relative elongation rate')
    if kw.get('dark'):
        m = m.replace('0.9     # k_d: light-triggered', '0.0     # k_d: dark')
    open(out, 'w').write(m)

if __name__ == '__main__':
    # 1. equilibrate: dark, growth frozen
    variant('hook.model', 'equil.model', nogrowth=True, dark=True)
    open('solver_eq.rk5','w').write('RK5Adaptive\n0 6\n0 2\n0.1 1e-4\n')
    subprocess.run([SIM, 'equil.model', 'hook.init', 'solver_eq.rk5',
                    '-init_output', 'hook_eq.init'],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    # 2. light + dark production runs from the equilibrated state
    variant('hook.model', 'hook_dark.model', dark=True)
    for model, out in (('hook.model','run_light.out'), ('hook_dark.model','run_dark.out')):
        subprocess.run([SIM, model, 'hook_eq.init', sys.argv[1] if len(sys.argv)>1 else 'solver0.rk5'],
                       stdout=open(out,'w'), stderr=subprocess.DEVNULL, check=True)
    print("done")
