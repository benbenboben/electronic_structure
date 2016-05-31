

import cython
from cython.parallel import prange
import numpy as np
cimport numpy as np


cpdef ao2mo_1body(int nb, double[:, :] hcore, double[:, :] cvec):
    cdef double[:, :] hcore_mo = np.zeros((nb, nb), dtype=float)
    cdef int p, q, mu, nu
    with nogil:
        for p in prange(nb):
            for q in range(nb):
                for mu in range(nb):
                    for nu in range(nb):
                        hcore_mo[p, q] += cvec[mu, p] * hcore[mu, nu] * cvec[nu, q]
    return hcore_mo