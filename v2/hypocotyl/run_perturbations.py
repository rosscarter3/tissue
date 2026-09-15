#!/usr/bin/env python3
"""Simulate the paper's pharmacological and genetic perturbations by changing
only the model parameter each treatment acts on, and nothing else.

  oryzalin   depolymerises microtubules -> no stress-guided fibre anisotropy,
             so the Eq. 3 redistribution is switched off (K_hill -> large,
             g(a) -> 0, fibre stays isotropic)
  isoxaben   inhibits cellulose synthesis -> less load-bearing fibre, so
             Y_fiber is reduced (dose-dependently)
  YUC6-OX    overproduces auxin -> the auxin gate is never released (k_d = 0)
  low light  weaker stimulus -> slower auxin depletion and weaker acidification
"""
import os, re, subprocess, sys
SIM = '/Users/ross/projects/tissue/v2/build/simulator'
BASE = open('hook.model').read()

def zero(m, tag, value='0.0'):
    return re.sub(r'^[0-9.eE+-]+(\s*#\s*@' + tag + r')', value + r'\1', m, flags=re.M)

TREATMENTS = {
    # name:        (edit function, description)
    'oryzalin':  lambda m: m.replace('0.7     # K_hill (paper: k = 0.7)',
                                     '50.0    # K_hill: MTs depolymerised, no guided anisotropy'),
    'isox100':   lambda m: m.replace('1350.0  # Y_fiber', '470.0   # Y_fiber (isoxaben 100 nM)'),
    'isox600':   lambda m: m.replace('1350.0  # Y_fiber', '200.0   # Y_fiber (isoxaben 600 nM)'),
    'yuc6':      lambda m: zero(m, 'AUXIN_DECAY'),
    'lowlight':  lambda m: zero(zero(m, 'AUXIN_DECAY', '0.15'), 'ACID_DRIVE', '0.7'),
}

if __name__ == '__main__':
    solver = sys.argv[1] if len(sys.argv) > 1 else 'solver.rk5'
    procs = []
    for name, edit in TREATMENTS.items():
        mf = f'pt_{name}.model'
        open(mf, 'w').write(edit(BASE))
        procs.append((name, subprocess.Popen(
            [SIM, mf, 'hook_eq.init', solver, '-centerTri_init'],
            stdout=open(f'pt_{name}.out', 'w'), stderr=subprocess.DEVNULL)))
    for name, p in procs:
        p.wait()
        print(f"{name}: exit {p.returncode}")
    print("perturbations done")
