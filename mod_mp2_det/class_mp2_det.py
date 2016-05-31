

import numpy as np
from copy import deepcopy
import itertools as it
from common.class_slater_condon import singlecomp_slater_condon
from common import helpers
from common import space_spin


class MP2Det(object):

    def __init__(self, nocc, nspace, ij, ab):
        self._nocc = nocc
        self._nspace = nspace
        self._nspin = 2 * self._nspace
        self._exc_level = len(ij)
        self._nocc_list = list(range(nocc))
        for i, a in zip(ij, ab):
            print(i, a)
            self._nocc_list[i] = a

    def __gen_list(self):
        low = max(self._exc_level - 2, 0)
        high = self._exc_level + 2
        for iexc in range(low, high + 1):
            if iexc == 0:
                tmpdet = list(range(self._nocc))
                if tmpdet != self._nocc_list:
                    yield tmpdet
            elif iexc == 1:
                for i in range(self._nocc):
                    for a in range(self._nocc, self._nspin):
                        tmpdet = list(range(self._nocc))
                        tmpdet[i] = a
                        if tmpdet != self._nocc_list:
                            yield tmpdet
            elif iexc == 2:
                for i, j in it.combinations(range(self._nocc), iexc):
                    for a, b in it.combinations(range(self._nocc, self._nspin), iexc):
                        tmpdet = list(range(self._nocc))
                        tmpdet[i] = a
                        tmpdet[j] = b
                        if tmpdet != self._nocc_list:
                            yield tmpdet
            elif iexc == 3:
                for i, j, k in it.combinations(range(self._nocc), iexc):
                    for a, b, c in it.combinations(range(self._nocc, self._nspin), iexc):
                        tmpdet = list(range(self._nocc))
                        tmpdet[i] = a
                        tmpdet[j] = b
                        tmpdet[k] = c
                        if tmpdet != self._nocc_list:
                            yield tmpdet

    def _calc_orbeng(self, hcore, vee):
        orbeng = np.zeros(hcore.shape[0], dtype=float)
        for i in range(self._nspin):
                orbeng[i] = hcore[i, i]
                for b in range(self._nocc):
                        orbeng[i] += vee[i, i, b, b] - vee[i, b, b, i]
        return orbeng

    @staticmethod
    def _calc_zeroth_eng(det, orbengs):
        eng = 0.0e0
        for i in det:
            eng += orbengs[i]
        return eng

    def driver(self, fhcore_mo, fvee_mo):
        hcore = space_spin.mo_space2spin_hcore(self._nspace, helpers.read_matrix_2d(fhcore_mo, self._nspace))
        vee = space_spin.mo_space2spin_vee(self._nspace, helpers.read_matrix_4d(fvee_mo, self._nspace))
        energy = 0.0e0
        orbeng = self._calc_orbeng(hcore, vee)
        E_i = self._calc_zeroth_eng(self._nocc_list, orbeng)
        for d in self.__gen_list():
            slater_condon = singlecomp_slater_condon(self._nocc_list, d)
            #mat_elem = slater_condon.calc_matrix_element(hcore, vee)
            two_body = slater_condon.calc_two_body(vee)
            hf_potential = slater_condon.calc_hf_potential(self._nocc, self._nspin, vee, self._nocc_list)
            mat_elem = two_body - hf_potential
            E_n = self._calc_zeroth_eng(d, orbeng)
            #print(type(two_body), type(hf_potential), type(E_i), type(E_n), type(mat_elem), type(mat_elem**2))
            if abs(E_i - E_n) < 1.0e-10: continue
            energy += (mat_elem**2) / (E_i - E_n)
        print('mp2: ', energy)
        slater_condon = singlecomp_slater_condon(self._nocc_list, self._nocc_list)
        two_body = slater_condon.calc_two_body(vee)
        hf_potential = slater_condon.calc_hf_potential(self._nocc, self._nspin, vee, self._nocc_list)
        first_order = two_body - hf_potential
        print('two body, hf_potential ',two_body, hf_potential)
        print('first order: ', first_order)
        print('zeroth order: ', E_i)
        print('total: ', energy + E_i + first_order)