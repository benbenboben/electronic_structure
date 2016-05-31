

import numpy as np
import scipy.linalg as la
from common import helpers


class SCF(object):

    def __init__(self, nocc, nspace):
        self._nocc = nocc
        self._nspace = nspace

    def driver(self, fhcore, fvee, fsmat, feigvec, feigval):
        smat = helpers.read_matrix_2d(fsmat, self._nspace)
        hcore_ao = helpers.read_matrix_2d(fhcore, self._nspace)
        vee_ao = helpers.read_matrix_4d(fvee, self._nspace)

        # calculate guess density matrix using hcore guess
        pmat = self._calc_pmat(hcore_ao, smat)

        #calculate norms of the eigenvectors
        norms = self._calc_norms(pmat)

        ncycle = 0
        while True:
            ncycle += 1
            fmat = self._calc_fmat(hcore_ao, vee_ao, pmat)
            pmat = self._calc_pmat(fmat, smat)
            norms2 = self._calc_norms(pmat)
            if max(abs(norms - norms2)) < 1.0e-10: break
            norms = norms2

        print("eng: {} ({} cycles)".format(self._calc_eng(hcore_ao, fmat, pmat), ncycle))

        eigvals, eigvecs = la.eigh(fmat, smat)

        helpers.print_matrix_2d(feigvec, eigvecs)
        helpers.print_vector_1d(feigval, eigvals)

    def _calc_transmat(self, smat):
        eigvals, eigvecs = la.eigh(smat)
        transmat = np.zeros((self._nspace, self._nspace), dtype=float)
        for i in range(transmat.shape[1]):
            transmat[:, i] = eigvecs[:, i] / eigvals[i]
        return transmat

    def _calc_pmat(self, fock, smat):
        eigvals, eigvecs = la.eigh(fock, smat)
        pmat = np.zeros((self._nspace, self._nspace), dtype=float)
        for mu in range(self._nspace):
            for nu in range(self._nspace):
                for a in range(self._nocc // 2):
                    pmat[mu, nu] += 2.0e0 * eigvecs[mu, a] * eigvecs[nu, a]
        return pmat

    def _calc_fmat(self, hcore_ao, vee_ao, pmat, scale_hcore=1.0e0, scale_vee=1.0e0):
        fmat = np.zeros((self._nspace, self._nspace), dtype=float)
        for mu in range(self._nspace):
            for nu in range(self._nspace):
                v = hcore_ao[mu, nu] * scale_hcore
                fmat[mu, nu] += v
                for lmda in range(self._nspace):
                    for sigma in range(self._nspace):
                        v = (vee_ao[mu, nu, sigma, lmda] - (0.5e0 * vee_ao[mu, lmda, sigma, nu])) * scale_vee
                        fmat[mu, nu] += pmat[lmda, sigma] * v
        return fmat

    def _calc_norms(self, pmat):
        norms = np.zeros(self._nspace, dtype=float)
        for i in range(self._nspace):
            norms[i] = la.norm(pmat[:, i])
        return norms

    @staticmethod
    def _unitary_transform(self, transmat, fmat):
        return np.dot(np.dot(transmat.T, fmat), transmat)

    def _calc_eng(self, hcore, fmat, pmat):
        eng = 0.0e0
        for mu in range(self._nspace):
            for nu in range(self._nspace):
                eng += pmat[mu, nu] * (hcore[mu, nu] + fmat[mu, nu])
        return eng * 0.5e0
