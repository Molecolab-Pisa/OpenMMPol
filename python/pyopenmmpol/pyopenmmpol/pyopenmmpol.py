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

_libopenmmpol.w_mmpol_init.argtypes = [ct.c_char_p]
_libopenmmpol.w_mmpol_init.restypes = []
def w_mmpol_init(infile_mmp):
    buf = ct.create_string_buffer(infile_mmp.encode())
    _libopenmmpol.w_mmpol_init(buf)

_libopenmmpol.do_mm.argtypes = []
_libopenmmpol.do_mm.restypes = []
def do_mm():
    _libopenmmpol.do_mm()

_libopenmmpol.do_qmmm.restypes = []
def do_qmmm(vqm, eqm):
    _eqm = np.ascontiguousarray(eqm)
    _vqm = np.ascontiguousarray(vqm)
    vqm_type = npct.ndpointer(dtype=ct.c_double,
                              ndim=2,
                              shape=(get_mm_atoms(),
                                     get_ld_cart()),
                                     flags='C_CONTIGUOUS')
    eqm_type = npct.ndpointer(dtype=ct.c_double,
                              ndim=3,
                              shape=(get_n_ipd(),
                                     get_pol_atoms(),
                                     3),
                                     flags='C_CONTIGUOUS')
    _libopenmmpol.do_qmmm.argtypes = [vqm_type, 
                                      eqm_type,
                                      INT_TYPE,
                                      INT_TYPE,
                                      INT_TYPE,
                                      INT_TYPE]
                                                  
    _libopenmmpol.do_qmmm(_vqm, 
                          _eqm,
                          get_ld_cart(),
                          get_mm_atoms(),
                          get_pol_atoms(),
                          get_n_ipd())

_libopenmmpol.restart.argtypes = []
_libopenmmpol.restart.restypes = []
def restart():
    _libopenmmpol.restart()

_libopenmmpol.get_energy.argtypes = [ct.POINTER(ct.c_double), 
                                     ct.POINTER(ct.c_double)]
_libopenmmpol.get_energy.restypes = []
def get_energy():
    EMM = ct.c_double(0.0)
    EPol = ct.c_double(0.0)
    _libopenmmpol.get_energy(ct.byref(EMM), 
                             ct.byref(EPol))
    return EMM.value, EPol.value

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

_libopenmmpol.is_amoeba.argtypes = []
_libopenmmpol.is_amoeba.restype = ct.c_bool
def is_amoeba():
     return bool(_libopenmmpol.is_amoeba())

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

class MMPol(object):
  def __init__(self, mmpol_file):
   
    if (get_mm_atoms() > 0):
      raise( RuntimeError('\n MMPol Module was already initialized! '))

    w_mmpol_init(mmpol_file)

    pass

  @property
  def is_amoeba(self):
    return is_amoeba()

  @property
  def n_ipd(self):
    return get_n_ipd()

  @property
  def ld_cart(self):
    return get_ld_cart()

  @property
  def mm_atoms(self):
    return get_mm_atoms()

  @property
  def pol_atoms(self):
    return get_pol_atoms()

  @property
  def cmm(self):
    return get_cmm()

  @property
  def cpol(self):
    return get_cpol()

  @property
  def q(self):
    return get_q()

  @property
  def ipd(self):
    return get_ipd()


  def do_mm(self):
    do_mm()

    
  def restart(self):
    restart()
    

  def do_qmmm(self,VQM,EQM):

    if VQM.ndim == 1:
      VQM = VQM[:,None] # ld_cart = 1
    if EQM.ndim == 2:
      EQM = EQM[None,:] # n_ipd = 1
    
    # Check shapes
    assert np.shape(VQM) == (self.mm_atoms,self.ld_cart)
    assert np.shape(EQM) == (self.n_ipd,self.pol_atoms,3)

    do_qmmm(VQM,EQM)

  def get_MM_energy(self):

    self.EMM, self.EMMPol = get_energy()

    print(" ---------------------------------------------")
    print("   MM     -- MM (q)    {:14.4f} kcal/mol".format(self.EMM*627.5094740630558))
    print("   MM     -- MM (dip)  {:14.4f} kcal/mol".format(self.EMMPol*627.5094740630558))
    print(" Total MM Energy       {:14.4f} kcal/mol".format((self.EMM+self.EMMPol)*627.5094740630558))

    return self.EMM,self.EMMPol

  def get_QMMM_energy(self,VQM,EQM):

    if VQM.ndim == 1:
      VQM = VQM[:,None] # ld_cart = 1
    if EQM.ndim == 2:
      EQM = EQM[None,:] # n_ipd = 1
    
    # Check shapes
    assert np.shape(VQM) == (self.mm_atoms,self.ld_cart)
    assert np.shape(EQM) == (self.n_ipd,self.pol_atoms,3)

    # Potential * multipoles
    self.e_static = np.sum(VQM*self.q)

    # Field * induced dipoles 
    self.e_polar  = -0.5 * np.sum(EQM*self.ipd)

    e = self.e_static + self.e_polar
    print(" ---------------------------------------------")
    print("   QM     -- MM (q)    {:14.4f} kcal/mol".format(self.e_static*627.5094740630558))
    print("   QM     -- MM (dip)  {:14.4f} kcal/mol".format(self.e_polar*627.5094740630558))
    print(" Total QM/MM Energy    {:14.4f} kcal/mol".format(e*627.5094740630558))


    return self.e_static,self.e_polar


#### For testing
# if __name__ == '__main__':
# 
#   # Initialize mmpol
#   mmp = MMPol('test01.mmp') 
# 
# 
#   # Get coordinates of MM part
#   cmm  = mmp.cmm
#   cpol = mmp.cpol
# 
#   # Do MM part
#   mmp.do_mm()
# 
#   # Fake a potential and a field
#   # (they should be calculated from QM density ...)
#   V = -cmm[:,0]*10
#   E = np.zeros(cpol.shape)
#   E[:,0] = 1
# 
# #  EQMMM = np.dot(V,mmp.q[:,0])
# #  print( EQMMM)
# 
#   print (V[:,None].T.shape)
# 
#   mmp.do_qmmm(V[:,None],E)
# 
#   print( mmpol.ipd.T )
# 
#   mmp.get_MM_energy()
#   mmp.get_QMMM_energy(V[:,None],E)
#   quit()

