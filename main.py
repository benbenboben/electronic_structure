

from mod_scf.class_scf import SCF
from mod_transform.transform import Transform
from mod_mp2.mp2 import MP2
from mod_mp2_det.class_mp2_det import MP2Det
from common.filenames import *


def main():
    nocc = 10
    nspace = 19

    scf = SCF(nocc, nspace)
    scf.driver(FILE_HCORE_AO, FILE_VEE_AO, FILE_SMAT_AO, FILE_EIGVEC, FILE_EIGVAL)

    transform = Transform(nocc, nspace)
    transform.driver(FILE_HCORE_AO, FILE_VEE_AO, FILE_EIGVEC, FILE_HCORE_MO, FILE_VEE_MO)

    #mp2 = MP2(nocc, nspace)
    #mp2.driver(FILE_HCORE_MO, FILE_VEE_MO)

    mp2det = MP2Det(nocc, nspace, [], [])
    mp2det.driver(FILE_HCORE_MO, FILE_VEE_MO)

    #mp2det2 = MP2Det(nocc, nspace, [9], [11])
    #mp2det2.driver(FILE_HCORE_MO, FILE_VEE_MO)

if __name__ == '__main__':
    main()
