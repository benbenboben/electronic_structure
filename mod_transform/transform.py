

import numpy as np
from common import helpers
from common import space_spin
import datetime


class Transform(object):

    def __init__(self, nocc, nspace):
        self._nocc = nocc
        self._nspace = nspace

    def driver(self, fhcore, fvee, feigvec, fhcore_mo, fvee_mo):
        hcore_ao = helpers.read_matrix_2d(fhcore, self._nspace)
        vee_ao = helpers.read_matrix_4d(fvee, self._nspace)
        eigvec = helpers.read_matrix_2d(feigvec, self._nspace, sym=False)
        hcore_mo_space = self._ao2mo_hcore(hcore_ao, eigvec)
        vee_mo_space = self._ao2mo_vee2(vee_ao, eigvec)
        #hcore_mo_spin = space_spin.mo_space2spin_hcore(self._nspace, hcore_mo_space)
        #vee_mo_spin = space_spin.mo_space2spin_vee(self._nspace, vee_mo_space)
        #self._calc_eng(hcore_mo_spin, vee_mo_spin)
        helpers.print_matrix_2d(fhcore_mo, hcore_mo_space)
        helpers.print_matrix_4d(fvee_mo, vee_mo_space)

    def _ao2mo_hcore(self, hcore, cvec):
        hcore_mo = np.dot(np.dot(cvec.T,hcore), cvec)
        return hcore_mo
        #hcore_mo = np.zeros((self._nspace, self._nspace), dtype=float)
        #for p in range(self._nspace):  ## mo idx
        #    for q in range(self._nspace):  ## mo idx
        #        for mu in range(self._nspace):  ## ao idx
        #            for nu in range(self._nspace):  ## ao idx
        #                hcore_mo[p, q] += (cvec[mu, p] * hcore[mu, nu] * cvec[nu, q])
        #                hcore_mo[q, p] = hcore_mo[p, q]
        #return hcore_mo
        #return ao2mo.ao2mo_1body(hcore.shape[0], hcore, cvec)

    def _ao2mo_vee2(self, vee, cvec):
        temp1 = np.einsum('ijkl,lz->ijkz', vee, cvec)
        temp2 = np.einsum('ijkl,kz->ijzl', temp1, cvec)
        temp3 = np.einsum('ijkl,jz->izkl', temp2, cvec)
        veemo = np.einsum('ijkl,iz->zjkl', temp3, cvec)
        return veemo

    def _ao2mo_vee(self, vee, cvec):
        nspace = self._nspace
        veemo = np.zeros((nspace, nspace, nspace, nspace), dtype=float)
        temp1 = np.zeros((nspace, nspace, nspace, nspace), dtype=float)
        temp2 = np.zeros((nspace, nspace, nspace, nspace), dtype=float)
        temp3 = np.zeros((nspace, nspace, nspace, nspace), dtype=float)
        for s in range(nspace):
            for mu in range(nspace):
                for nu in range(nspace):
                    for rho in range(nspace):
                        for sgma in range(nspace):
                            temp1[mu, nu, rho, s] += vee[mu, nu, rho, sgma] * cvec[sgma, s]
        for r in range(nspace):
            for s in range(nspace):
                for mu in range(nspace):
                    for nu in range(nspace):
                        for rho in range(nspace):
                            temp2[mu, nu, r, s] += temp1[mu, nu, rho, s] * cvec[rho, r]
        for q in range(nspace):
            for r in range(nspace):
                for s in range(nspace):
                    for mu in range(nspace):
                        for nu in range(nspace):
                            temp3[mu, q, r, s] += temp2[mu, nu, r, s] * cvec[nu, q]
        for p in range(nspace):
            for q in range(nspace):
                for r in range(nspace):
                    for s in range(nspace):
                        for mu in range(nspace):
                            veemo[p, q, r, s] += temp3[mu, q, r, s] * cvec[mu, p]
        return veemo


    def _calc_eng(self, hcore_mo_spin, vee_mo_spin):
        eng1 = 0.0
        eng2 = 0.0
        for i in range(self._nocc):
            eng1 += hcore_mo_spin[i, i]
            for j in range(self._nocc):
                eng2 += 0.5 * (vee_mo_spin[i, i, j, j] - vee_mo_spin[i, j, j, i])
        print(eng1 + eng2)