

import numpy as np


def upper_pack(i, j):
    """
    i: type(int): index
    j: type(int): index
    returns: type(int): UP index of i, j
    """
    min_idx = min(i, j)
    max_idx = max(i, j)
    return (max_idx * (max_idx + 1) // 2) + min_idx


def spin_of(i):
    """
    i: type(int): index
    returns: type(int): returns 0 or 1 (alpha or beta spin)
    """
    if i % 2 == 0:
        return 0
    else:
        return 1


def spin_to_space(i):
    """
    i: type(int): index
    returns: type(int): SPATIAL index for SPIN index i
    """
    if i % 2 == 0:
        return i // 2
    else:
        return ((i + 1) // 2) - 1


def read_matrix_2d(fname, size, idx='f', sym=True):
    """
    fname: type(str): name of matrix file
    returns: numpy ndarray of matrix
    """
    with open(fname) as f:
        lines = f.readlines()
    mat = np.zeros((size, size), dtype=float)
    for i, l in enumerate(lines):
        if i == 0: continue
        l = l.strip().split()
        i = int(l[0])
        j = int(l[1])
        if i > size or j > size:
            continue
        if idx == 'f':
            i -= 1
            j -= 1
        v = float(l[2])
        mat[i, j] = v
        if sym:
            mat[j, i] = v
    return mat

def read_matrix_4d(fname, size, idx='f'):
    """
    fname: type(str): name of matrix file
    returns: numpy ndarray of matrix
    """
    with open(fname) as f:
        lines = f.readlines()
    mat = np.zeros((size, size, size, size), dtype=float)
    for i, l in enumerate(lines):
        if i == 0: continue
        l = l.strip().split()
        mu = int(l[0])
        nu = int(l[1])
        lmda = int(l[2])
        sgma = int(l[3])
        if mu > size or nu > size or lmda > size or sgma > size:
            continue
        if idx == 'f':
            mu -= 1
            nu -= 1
            lmda -= 1
            sgma -= 1
        val = float(l[4])
        mat[mu, nu, lmda, sgma] = \
        mat[nu, mu, lmda, sgma] = \
        mat[mu, nu, sgma, lmda] = \
        mat[nu, mu, sgma, lmda] = \
        mat[lmda, sgma, mu, nu] = \
        mat[sgma, lmda, mu, nu] = \
        mat[lmda, sgma, nu, mu] = \
        mat[sgma, lmda, nu, mu] = \
        val
    return mat


def print_matrix_2d(fname, mat, idx='f'):
    with open(fname, 'w') as f:
        f.write('{}    {}'.format(mat.shape[0], mat.shape[1]))
        f.write('\n')
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if idx == 'f':
                    f.write('{}    {}    {}'.format(i + 1, j + 1, mat[i, j]))
                else:
                    f.write('{}    {}    {}'.format(i, j, mat[i, j]))
                f.write('\n')
    f.close()


def print_matrix_4d(fname, mat, idx='f'):
    with open(fname, 'w') as f:
        f.write('{}    {}'.format(mat.shape[0], mat.shape[2]))
        f.write('\n')
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                for k in range(mat.shape[2]):
                    for l in range(mat.shape[3]):
                        if idx == 'f':
                            f.write('{}    {}    {}    {}    {}'.format(i + 1, j + 1, k + 1, l + 1, mat[i, j, k, l]))
                        else:
                            f.write('{}    {}    {}    {}    {}'.format(i, j, k ,l , mat[i, j, k, l]))
                        f.write('\n')
    f.close()


def print_vector_1d(fname, v, idx='f'):
    with open(fname, 'w') as f:
        f.write('{}'.format(v.shape[0]))
        f.write('\n')
        for i in range(v.shape[0]):
            if idx == 'f':
                f.write('{}    {}'.format(i + 1, v[i]))
            else:
                f.write('\n')
    f.close()


def read_vector_1d(fname, size, idx='f'):
    with open(fname) as f:
        lines = f.readlines()
    mat = np.zeros(size, dtype=float)
    for i, l in enumerate(lines):
        if i == 0: continue
        l = l.strip().split()
        i = int(l[0])
        if i > size:
            continue
        if idx == 'f':
            i -= 1
        v = float(l[1])
        mat[i] = v
    return mat

