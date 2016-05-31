

import numpy as np


def mo_space2spin_hcore(nspace, hcore):
    nspin = 2 * nspace
    spinmat = make_spinmat(nspace)
    hcore_spin = np.zeros((nspin, nspin), dtype=float)
    for ispc in range(nspace):
        for ispin in range(2):
            imo = spinmat[ispc, ispin]
            for jspc in range(nspace):
                for jspin in range(2):
                    jmo = spinmat[jspc, jspin]
                    if ispin != jspin:
                        continue
                    hcore_spin[imo, jmo] = hcore[ispc, jspc]
    return hcore_spin


def mo_space2spin_vee(nspace, vee):
    nspin = 2 * nspace
    spinmat = make_spinmat(nspace)
    vee_spin = np.zeros((nspin, nspin, nspin, nspin), dtype=float)
    for ispc in range(nspace):
        for ispin in range(2):
            imo = spinmat[ispc, ispin]
            for jspc in range(nspace):
                for jspin in range(2):
                    jmo = spinmat[jspc, jspin]
                    if ispin != jspin:
                        continue
                    for kspc in range(nspace):
                        for kspin in range(2):
                            kmo = spinmat[kspc, kspin]
                            for lspc in range(nspace):
                                for lspin in range(2):
                                    lmo = spinmat[lspc, lspin]
                                    if kspin != lspin:
                                        continue
                                    vee_spin[imo, jmo, kmo, lmo] = vee[ispc, jspc, kspc, lspc]
    return vee_spin


def make_spinmat(nspace):
    spinmat = np.zeros((nspace, 2), dtype=int)
    imo = -1
    for ispc in range(nspace):
        for ispn in range(2):
            imo += 1
            spinmat[ispc, ispn] = imo
    return spinmat