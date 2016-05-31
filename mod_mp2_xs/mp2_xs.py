

import numpy as np
from common import space_spin
from common import helpers


class MP2XS(object):

    def __init__(self, nocc, nspace, ij, ab):
        self._nocc = nocc
        self._nspace = nspace
        self._nspin = self._nspace * 2
        self._nocc_list = list(range(self._nocc))
        self._virt_list = list(range(self._nspin - self._nocc))
        for i, a in zip(ij, ab):
            print(i, a)
            self._nocc_list[i] = a
            del self._virt_list[self._virt_list.index(a)]

    def driver(self, fhcore_mo, fvee_mo):
        hcore = space_spin.mo_space2spin_hcore(self._nspace, helpers.read_matrix_2d(fhcore_mo, self._nspace))
        vee = space_spin.mo_space2spin_vee(self._nspace, helpers.read_matrix_4d(fvee_mo, self._nspace))
        orb_eng = self._calc_orbeng(hcore, vee)
        #tvec = self._calc_tvec(vee, orb_eng)
        print(self._calc_eng2(orb_eng, vee))

    def _construct_fock(self, hcore, vee):
        fmat = np.zeros(hcore.shape, dtype=float)
        for i in range(self._nspin):
            for j in range(self._nspin):
                fmat[i, j] = hcore[i, j]
                for a in range(self._nspin):
                        fmat[i, j] += vee[i, j, a, a] - vee[i, a, j, a]
        return fmat

    def _calc_orbeng(self, hcore, vee):
        fmat = np.zeros(hcore.shape[0], dtype=float)
        for i in range(self._nspin):
                fmat[i] = hcore[i, i]
                for b in range(self._nocc):
                        fmat[i] += vee[i, i, b, b] - vee[i, b, b, i]
        return fmat

    def _calc_tvec(self, vee, orb_eng):
        tvec = np.zeros((self._nocc, self._nocc, self._nspin, self._nspin), dtype=float)
        for i in range(self._nocc):
            for j in range(self._nocc):
                if i == j: continue
                for a in range(self._nocc, self._nspin):
                    for b in range(self._nocc, self._nspin):
                        if a == b: continue
                        eijab = orb_eng[i] + orb_eng[j] - orb_eng[a] - orb_eng[b]
                        intval = vee[i, a, j, b] - vee[i, b, j, a]
                        tvec[i, j, a, b] = intval / eijab
        return tvec

    def _calc_eng(self, tvec, vee):
        #return np.einsum('ijab,ijab', vee[:self._nocc, :self._nocc, :, :], tvec)
        #return np.tensordot(tvec, vee[:self._nocc, :self._nocc, :, :])
        sum = 0.e0
        for i in range(self._nocc):
            for j in range(self._nocc):
                for a in range(self._nocc, self._nspin):
                    for b in range(self._nocc, self._nspin):
                        sum += tvec[i, j, a, b] * (vee[i, a, j, b] - vee[i, b, j, a])
        return 0.25e0 * sum

    def _calc_eng2(self, orb_eng, vee):
        sum = 0.e0
        for i in self._nocc_list:
            for j in self._nocc_list:
                if i == j: continue
                for a in self._virt_list:
                    for b in self._virt_list:
                        if a == b: continue
                        numerator = (vee[i, a, j, b] - vee[i, b, j, a])**2
                        denominator = orb_eng[i] + orb_eng[j] - orb_eng[a] - orb_eng[b]
                        if abs(denominator) < 1.0e-10:
                            print('{}    {}    {}    {}    {}    {}'.format(i, j, a, b, numerator, denominator))
                            continue
                        try:
                            sum += numerator / denominator
                        except:
                            pass
        return 0.25e0 * sum