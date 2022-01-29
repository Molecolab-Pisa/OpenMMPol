"""
Python MMPol library

High-level python interface to the FORTRAN/f2py code

"""

import pymmpol
from pymmpol import mmpol # fortran module
import numpy as np


class MMPol(object):


  def __init__(self,mmpol_file):
    
    if (mmpol.mm_atoms > 0):
      raise( RuntimeError('\n MMPol Module was already initialized! '))

    pymmpol.w_mmpol_init(mmpol_file)

    pass

  @property
  def is_amoeba(self):
    return bool(mmpol.amoeba)

  @property
  def n_ipd(self):
    return mmpol.n_ipd

  @property
  def ld_cart(self):
    return mmpol.ld_cart

  @property
  def mm_atoms(self):
    return mmpol.mm_atoms

  @property
  def pol_atoms(self):
    return mmpol.pol_atoms

  @property
  def cmm(self): 
    return mmpol.cmm.T

  @property
  def cpol(self):
    return mmpol.cpol.T

  @property
  def q(self):
    return mmpol.q.T

  @property
  def ipd(self):
    return mmpol.ipd.T


  def do_mm(self):
    pymmpol.do_mm()


  def do_qmmm(self,VQM,EQM):

    if VQM.ndim == 1:
      VQM = VQM[:,None] # ld_cart = 1
    if EQM.ndim == 2:
      EQM = EQM[None,:] # n_ipd = 1
    
    # Check shapes
    assert np.shape(VQM) == (self.mm_atoms,self.ld_cart)
    assert np.shape(EQM) == (self.n_ipd,self.pol_atoms,3)

    pymmpol.do_qmmm(VQM.T,EQM.T)

  def get_MM_energy(self):

    self.EMM,self.EMMPol = pymmpol.get_energy()

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
if __name__ == '__main__':

  # Initialize mmpol
  mmp = MMPol('test01.mmp') 


  # Get coordinates of MM part
  cmm  = mmp.cmm
  cpol = mmp.cpol

  # Do MM part
  mmp.do_mm()

  # Fake a potential and a field
  # (they should be calculated from QM density ...)
  V = -cmm[:,0]*10
  E = np.zeros(cpol.shape)
  E[:,0] = 1

#  EQMMM = np.dot(V,mmp.q[:,0])
#  print( EQMMM)

  print (V[:,None].T.shape)

  mmp.do_qmmm(V[:,None],E)

  print( mmpol.ipd.T )

  mmp.get_MM_energy()
  mmp.get_QMMM_energy(V[:,None],E)
  quit()

