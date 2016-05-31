

from mod_scf.class_scf import SCF


def main():
    file_hcore_ao = 'hcore1.fmt'
    file_vee_ao = 'vee1.fmt'
    file_smat_ao = 'smat1.fmt'

    int = 2
    nspace = 2
    scf = SCF(int, nspace)
    scf.driver(file_hcore_ao, file_vee_ao, file_smat_ao)


if __name__ == '__main__':
    main()