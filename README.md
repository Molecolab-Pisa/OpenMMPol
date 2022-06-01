# Open-MMPol

## Description

Open-MMPol is a ibrary to interface quantum chemical software with atomistic polarizable embedding. 
It is written and mantained by the MoLECoLab (Modeling Light & Environment in Complex Systems) 
research group at the University of Pisa. 

## Installation

The package is mantained and released with git. To get the last version clone the repository:

``git clone git@molimen1.dcci.unipi.it:molecolab/open-mmpol.git``

To install the package, you just need a compiler that supports standard Fortran 2003, cmake, 
lapack libraries and optionally hdf5 library; to install the requirements on OpenSuse Leap just 
use the following command:

``zypper --non-interactive in gcc gcc-c++ gcc-fortran make cmake python lapack-devel liblapack3 hdf5 hdf5-devel zlib-devel``

To build the package, just use cmake in the standard way; from the root directory of the git repository:

``$ mkdir build``
  
``$ cd build``
  
``$ cmake -DCMAKE_BUILD_TYPE=DEBUG .. # You can use DEBUG or RELEASE to control the build flags``
  
``$ make``




