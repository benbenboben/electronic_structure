

from common import helpers
import copy


class singlecomp_slater_condon(object):

    def __init__(self, bra, ket):
        self._bra = copy.deepcopy(bra)
        self._ket = copy.deepcopy(ket)
        self.__max_coincidence()

    def calc_ovlap(self):
        for i in self._bra:
            if i not in self._ket:
                return 0.0e0
        else:
            return self._sign
        #if np.array_equal(self._bra, self._ket):
        #    return self._sign
        #else:
        #    return 0.0e0

    def calc_matrix_element(self, hcore, vee):
        ndiff, unaligned_indices = self.__calc_wrong_indices(self._bra, self._ket)
        one_body = self._calc_one_body(ndiff, unaligned_indices, hcore)
        two_body = self._calc_two_body(ndiff, unaligned_indices, vee)
        return one_body + two_body

    def calc_one_body(self, hcore):
        ndiff, unaligned_indices = self.__calc_wrong_indices(self._bra, self._ket)
        return self._calc_one_body(ndiff, unaligned_indices, hcore)

    def calc_two_body(self, vee):
        ndiff, unaligned_indices = self.__calc_wrong_indices(self._bra, self._ket)
        return self._calc_two_body(ndiff, unaligned_indices, vee)

    def _calc_one_body(self, ndiff, unaligned_indices, hcore):
        val = 0.0e0
        if ndiff == 0:
            for i in range(len(self._bra)):
                val += hcore[self._bra[i], self._ket[i]]
        elif ndiff == 1:
            idx = unaligned_indices[0]
            spin_scale = self.__scale_for_spin_one_body(self._bra[idx], self._ket[idx])
            val += (hcore[self._bra[idx], self._ket[idx]] * spin_scale)
        return val * self._sign

    def calc_hf_potential(self, nocc, nspin, vee, occlist):
        ndiff, unaligned_indices = self.__calc_wrong_indices(self._bra, self._ket)
        return self._calc_hf_potential(nocc, nspin, ndiff, unaligned_indices, vee, occlist)

    def _calc_hf_potential(self, nocc, nspin, ndiff, unaligned_indices, vee, occlist):
        val = 0.0e0
        if ndiff == 0:
            for i in range(len(self._bra)):
                #for a in range(nocc):
                for a in occlist:
                    scale_coul = scale_exch = 1.0e0
                    #scale_coul = self.__scale_for_spin_two_body(self._bra[i], self._ket[i], self._bra[a], self._ket[a])
                    #scale_exch = self.__scale_for_spin_two_body(self._bra[i], self._ket[a], self._bra[a], self._ket[i])
                    #val += (vee[self._bra[i], self._ket[i], self._bra[a], self._ket[a]] * scale_coul) \
                    #     - (vee[self._bra[i], self._ket[a], self._bra[a], self._ket[i]] * scale_exch)
                    val += (vee[self._bra[i], self._ket[i], a, a] * scale_coul) \
                         - (vee[self._bra[i], a, a, self._ket[i]] * scale_exch)
        elif ndiff == 1:
            idx = unaligned_indices[0]
            #for a in range(nocc):
            for a in occlist:
                scale_coul = scale_exch = 1.0e0
                #scale_coul = self.__scale_for_spin_two_body(self._bra[idx], self._ket[idx], self._bra[a], self._ket[a])
                #scale_exch = self.__scale_for_spin_two_body(self._bra[idx], self._ket[a], self._bra[a], self._ket[idx])
                #val += (vee[self._bra[idx], self._ket[idx], self._bra[a], self._ket[a]] * scale_coul) \
                #     - (vee[self._bra[idx], self._ket[a], self._bra[a], self._ket[idx]] * scale_exch)
                val += (vee[self._bra[idx], self._ket[idx], a, a] * scale_coul) \
                     - (vee[self._bra[idx], a, a, self._ket[idx]] * scale_exch)
        return val * self._sign

    def _calc_two_body(self, ndiff, unaligned_indices, vee):
        val = 0.0e0
        if ndiff == 0:
            for i in range(len(self._bra)):
                for j in range(i + 1, len(self._bra)):
                    scale_coul = self.__scale_for_spin_two_body(self._bra[i], self._ket[i], self._bra[j], self._ket[j])
                    scale_exch = self.__scale_for_spin_two_body(self._bra[i], self._ket[j], self._bra[j], self._ket[i])
                    val += (vee[self._bra[i], self._ket[i], self._bra[j], self._ket[j]] * scale_coul) \
                         - (vee[self._bra[i], self._ket[j], self._bra[j], self._ket[i]] * scale_exch)
        if ndiff == 1:
            idx = unaligned_indices[0]
            for i in range(len(self._bra)):
                scale_coul = self.__scale_for_spin_two_body(self._bra[idx], self._ket[idx], self._bra[i], self._ket[i])
                scale_exch = self.__scale_for_spin_two_body(self._bra[idx], self._ket[i], self._bra[i], self._ket[idx])
                val += (vee[self._bra[idx], self._ket[idx], self._bra[i], self._ket[i]] * scale_coul) \
                     - (vee[self._bra[idx], self._ket[i], self._bra[i], self._ket[idx]] * scale_exch)
        if ndiff == 2:
            idx1, idx2 = unaligned_indices[0], unaligned_indices[1]
            scale_coul = self.__scale_for_spin_two_body(self._bra[idx1], self._ket[idx1], self._bra[idx2], self._ket[idx2])
            scale_exch = self.__scale_for_spin_two_body(self._bra[idx1], self._ket[idx2], self._bra[idx2], self._ket[idx1])
            val += (vee[self._bra[idx1], self._ket[idx1], self._bra[idx2], self._ket[idx2]] * scale_coul) \
                 - (vee[self._bra[idx1], self._ket[idx2], self._bra[idx2], self._ket[idx1]] * scale_exch)
        return val * self._sign

    def __calc_ndiff(self):
        ndiff = 0
        for idx in self._bra:
            if idx not in self._ket:
                ndiff += 1
        return ndiff

    def __scale_for_spin_two_body(self, i, j, a, b):
        return self.__scale_for_spin_one_body(i, j) * self.__scale_for_spin_one_body(a, b)

    def __scale_for_spin_one_body(self, i, a):
        if helpers.spin_of(i) == helpers.spin_of(a):
            return 1.0e0
        else:
            return 0.0e0

    def __calc_wrong_indices(self, bra, ket):
        wrongs = 0
        indices = []
        for i in range(len(bra)):
            if bra[i] != ket[i]:
                wrongs += 1
                indices.append(i)
        return wrongs, indices

    def __max_coincidence(self):
        n = len(self._bra)
        sign = 1
        min_wrongs, _ = self.__calc_wrong_indices(self._bra, self._ket)
        for k in range(n * n):
            nswap = 0
            for i in range(n):
                for j in range(i + 1, n):
                    tmpket = self.__dpcpy(self._ket[:])
                    #tmpket = copy.deepcopy(self._ket)
                    tmp = tmpket[i]
                    tmpket[i] = tmpket[j]
                    tmpket[j] = tmp
                    curr_wrongs, _ = self.__calc_wrong_indices(self._bra, tmpket)
                    if curr_wrongs < min_wrongs:
                        min_wrongs = curr_wrongs
                        sign *= -1
                        self._ket = self.__dpcpy(tmpket[:])
                        #self._ket = copy.deepcopy(tmpket)
                        nswap += 1
            if nswap == 0:
                self._sign = sign
                break
    def __dpcpy(self, a):
        b = []
        for i in range(len(a)):
            b.append(a[i])
        return b
