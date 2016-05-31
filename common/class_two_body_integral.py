

import numpy as np

from common import helpers


class two_body_integral(object):

    def __init__(self, nocc, nspace, fname, scale=1.0e0, idx_lang='f'):
        """
        nocc: type(int): occupancy of system
        nspace:  type(int): number of virtuals of system (spatial orbitals)
        fname: tpye(str): file where integrals can be read from.
                          indices are expected to be spatial, chemist notation
        """
        assert type(nocc) is int
        assert type(nspace) is int
        assert type(fname) is str

        self._npack_half = helpers.upper_pack(nspace, nspace)
        self._npack_full = helpers.upper_pack(self._npack_half, self._npack_half)
        self._int_vals = np.zeros(self._npack_full, dtype=float)

        with open(fname, 'r') as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if i == 0: continue
            line = line.strip().split()
            p = int(line[0])
            q = int(line[1])
            s = int(line[2])
            r = int(line[3])
            if p > nspace or q > nspace or s > nspace or r > nspace:
                continue
            if idx_lang == 'f':
                p -= 1
                q -= 1
                s -= 1
                r -= 1
            v = float(line[4])
            pq = helpers.upper_pack(p, q)
            sr = helpers.upper_pack(s, r)
            pqsr = helpers.upper_pack(pq, sr)
            self._int_vals[pqsr] = v * scale

    def __getitem__(self, pqsr_tup):
        """
        retrieve integral values given SPIN occupancy pqsr
        pqsr_tup: type(tuple of int): SPIN index of particles
        returns: type(float): value of two body integral with spin integration
        """
        assert len(pqsr_tup) == 4
        p, q, s, r = pqsr_tup
        spin_int = helpers.spin_of(p) * helpers.spin_of(q) * helpers.spin_of(s) * helpers.spin_of(r)
        p_spatial = helpers.spin_to_space(p)
        q_spatial = helpers.spin_to_space(q)
        s_spatial = helpers.spin_to_space(s)
        r_spatial = helpers.spin_to_space(r)
        return self._int_vals[helpers.upper_pack(helpers.upper_pack(p_spatial, q_spatial),
                                                 helpers.upper_pack(s_spatial, r_spatial))] * spin_int
