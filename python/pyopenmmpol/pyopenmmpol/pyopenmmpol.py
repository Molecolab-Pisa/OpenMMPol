import numpy as np
import numpy.ctypeslib as npct
import ctypes as ct

try:
    _libopenmmpol = ct.CDLL('libopenmmpol.so')
except OSError:
    print("Cannot find libopenmmpol.so. Check to have set correctly $LD_LIBRARY_PATH")
    raise ImportError

if _libopenmmpol.__use_8bytes_int():
    INT_TYPE = ct.c_int64
else:
    INT_TYPE = ct.c_int32

_libopenmmpol.mmpol_init_mmp.argtypes = [ct.c_char_p]
_libopenmmpol.mmpol_init_mmp.restypes = []
def ommp_init_mmp(infile_mmp):
    buf = ct.create_string_buffer(infile_mmp.encode())
    _libopenmmpol.mmpol_init_mmp(buf)

_libopenmmpol.mmpol_init_xyz.argtypes = [ct.c_char_p]
_libopenmmpol.mmpol_init_xyz.restypes = []
def ommp_init_xyz(infile_prm, infile_xyz):
    bufprm = ct.create_string_buffer(infile_prm.encode())
    bufxyz = ct.create_string_buffer(infile_xyz.encode())
    _libopenmmpol.mmpol_init_xyz(bufxyz, bufprm)

_libopenmmpol.ommp_terminate.argtypes = []
_libopenmmpol.ommp_terminate.restypes = []
def ommp_terminate():
    _libopenmmpol.ommp_terminate()

# ------

def set_external_field(eqm, solver):
    _eqm = np.ascontiguousarray(eqm)
    eqm_type = npct.ndpointer(dtype=ct.c_double,
                              ndim=2,
                              shape=(get_pol_atoms(),
                                     3),
                              flags='C_CONTIGUOUS')

    _libopenmmpol.set_external_field.restypes = []
    _libopenmmpol.set_external_field.argtypes = [eqm_type, INT_TYPE]
    _libopenmmpol.set_external_field(_eqm, ???)


_libopenmmpol.get_polelec_energy.argtypes = [ct.POINTER(ct.c_double)]
_libopenmmpol.get_polelec_energy.restypes = []
def get_polelec_energy():
    EPol = ct.c_double(0.0)
    _libopenmmpol.get_polelec_energy(ct.byref(EPol))
    return EPol.value

_libopenmmpol.get_fixedelec_energy.argtypes = [ct.POINTER(ct.c_double)]
_libopenmmpol.get_fixedelec_energy.restypes = []
def get_fixedelec_energy():
    EMM = ct.c_double(0.0)
    _libopenmmpol.get_fixedelec_energy(ct.byref(EMM))
    return EMM.value

# -------

_libopenmmpol.get_vdw_energy.argtypes = [ct.POINTER(ct.c_double)]
_libopenmmpol.get_vdw_energy.restype = []
def get_vdw_energy():
    ev = ct.c_double(0.0)
    _libopenmmpol.get_vdw_energy(ct.byref(ev))
    return ev

_libopenmmpol.get_bond_energy.argtypes = [ct.POINTER(ct.c_double)]
_libopenmmpol.get_bond_energy.restype = []
def get_bond_energy():
    eb = ct.c_double(0.0)
    _libopenmmpol.get_bond_energy(ct.byref(eb))
    return eb

_libopenmmpol.get_angle_energy.argtypes = [ct.POINTER(ct.c_double)]
_libopenmmpol.get_angle_energy.restype = []
def get_angle_energy():
    ea = ct.c_double(0.0)
    _libopenmmpol.get_angle_energy(ct.byref(ea))
    return ea

_libopenmmpol.get_angtor_energy.argtypes = [ct.POINTER(ct.c_double)]
_libopenmmpol.get_angtor_energy.restype = []
def get_angtor_energy():
    eat = ct.c_double(0.0)
    _libopenmmpol.get_angtor_energy(ct.byref(eat))
    return eat

_libopenmmpol.get_strtor_energy.argtypes = [ct.POINTER(ct.c_double)]
_libopenmmpol.get_strtor_energy.restype = []
def get_strtor_energy():
    ebt = ct.c_double(0.0)
    _libopenmmpol.get_strtor_energy(ct.byref(ebt))
    return ebt

_libopenmmpol.get_strbnd_energy.argtypes = [ct.POINTER(ct.c_double)]
_libopenmmpol.get_strbnd_energy.restype = []
def get_strbnd_energy():
    eba = ct.c_double(0.0)
    _libopenmmpol.get_strbnd_energy(ct.byref(eba))
    return eba

_libopenmmpol.get_opb_energy.argtypes = [ct.POINTER(ct.c_double)]
_libopenmmpol.get_opb_energy.restype = []
def get_opb_energy():
    eopb = ct.c_double(0.0)
    _libopenmmpol.get_opb_energy(ct.byref(eopb))
    return eopb

_libopenmmpol.get_pitors_energy.argtypes = [ct.POINTER(ct.c_double)]
_libopenmmpol.get_pitors_energy.restype = []
def get_pitors_energy():
    ept = ct.c_double(0.0)
    _libopenmmpol.get_pitors_energy(ct.byref(ept))
    return ept

_libopenmmpol.get_torsion_energy.argtypes = [ct.POINTER(ct.c_double)]
_libopenmmpol.get_torsion_energy.restype = []
def get_torsion_energy():
    et = ct.c_double(0.0)
    _libopenmmpol.get_torsion_energy(ct.byref(et))
    return et

_libopenmmpol.get_tortor_energy.argtypes = [ct.POINTER(ct.c_double)]
_libopenmmpol.get_tortor_energy.restype = []
def get_tortor_energy():
    ett = ct.c_double(0.0)
    _libopenmmpol.get_tortor_energy(ct.byref(et))
    return ett

_libopenmmpol.get_urey_energy.argtypes = [ct.POINTER(ct.c_double)]
_libopenmmpol.get_urey_energy.restype = []
def get_urey_energy():
    eub = ct.c_double(0.0)
    _libopenmmpol.get_urey_energy(ct.byref(eub))
    return eub

# -------

_libopenmmpol.get_n_ipd.argtypes = []
_libopenmmpol.get_n_ipd.restype = INT_TYPE
def get_n_ipd():
    return int(_libopenmmpol.get_n_ipd())

_libopenmmpol.get_ld_cart.argtypes = []
_libopenmmpol.get_ld_cart.restype = INT_TYPE
def get_ld_cart():
    return int(_libopenmmpol.get_ld_cart())

_libopenmmpol.get_mm_atoms.argtypes = []
_libopenmmpol.get_mm_atoms.restype = INT_TYPE
def get_mm_atoms():
    return int(_libopenmmpol.get_mm_atoms())

_libopenmmpol.get_pol_atoms.argtypes = []
_libopenmmpol.get_pol_atoms.restype = INT_TYPE
def get_pol_atoms():
    return int(_libopenmmpol.get_pol_atoms())

_libopenmmpol.get_cmm.argtypes = []
def get_cmm():
    rtype = npct.ndpointer(dtype=ct.c_double,
                           ndim=2,
                           shape=(get_mm_atoms(),3),
                           flags=('F_CONTIGUOUS'))
    _libopenmmpol.get_cmm.restype = rtype
    cmm = _libopenmmpol.get_cmm()
    return cmm

_libopenmmpol.get_cpol.argtypes = []
def get_cpol():
    rtype = npct.ndpointer(dtype=ct.c_double,
                           ndim=2,
                           shape=(get_pol_atoms(),3),
                           flags=('F_CONTIGUOUS'))
    _libopenmmpol.get_cpol.restype = rtype
    cpol = _libopenmmpol.get_cpol()
    return cpol

_libopenmmpol.get_q.argtypes = []
def get_q():
    rtype = npct.ndpointer(dtype=ct.c_double,
                           ndim=2,
                           shape=(get_mm_atoms(), get_ld_cart()),
                           flags=('F_CONTIGUOUS'))
    _libopenmmpol.get_q.restype = rtype
    q = _libopenmmpol.get_q()
    return q

_libopenmmpol.get_ipd.argtypes = []
def get_ipd():
    rtype = npct.ndpointer(dtype=ct.c_double,
                           ndim=3,
                           shape=(get_n_ipd(), get_pol_atoms(), 3),
                           flags=('F_CONTIGUOUS'))
    _libopenmmpol.get_ipd.restype = rtype
    ipd = _libopenmmpol.get_ipd()
    return ipd

# get_polar_mm

_libopenmmpol.is_amoeba.argtypes = []
_libopenmmpol.is_amoeba.restype = ct.c_bool
def is_amoeba():
     return bool(_libopenmmpol.is_amoeba())

