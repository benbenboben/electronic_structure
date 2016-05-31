

import numpy as np

from common import helpers


class one_body_integral(object):

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

        self._npack = helpers.upper_pack(nspace, nspace)
        self._int_vals = np.zeros(self._npack, dtype=float)

        with open(fname, 'r') as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if i == 0: continue
            line = line.strip().split()
            p = int(line[0])
            q = int(line[1])
            if p > nspace or q > nspace:
                continue
            if idx_lang == 'f':
                p -= 1
                q -= 1
            if len(line) == 4:
                v = float(line[4])
            else:
                v = float(line[2])
            pq = helpers.upper_pack(p, q)
            self._int_vals[pq] = v * scale

    def __getitem__(self, pq_tup):
        """
        retrieve integral values given SPIN occupancy p and q
        pq_tup: type(tuple of int): SPIN index of particle p, q
        returns: type(float): value of one body integral with spin integration
        """
        assert len(pq_tup) == 2
        p, q = pq_tup
        spin_int = helpers.spin_of(p) * helpers.spin_of(q)
        p_spatial = helpers.spin_to_space(p)
        q_spatial = helpers.spin_to_space(q)
        return self._int_vals[helpers.upper_pack(p_spatial, q_spatial)] * spin_int
