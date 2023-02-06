#ifndef __F_CART_COMPONENTS__
#define __F_CART_COMPONENTS__

#define _amoeba_D_ 1
!! Index of direct (D) field and dipoles in AMOEBA FF (true dipoles)
#define _amoeba_P_ 2
!! Index of polarization (P) field and dipoles in AMOEBA FF (Lagrange multiplier)

#define _x_ 1
#define _y_ 2
#define _z_ 3

#define _xx_ 1
#define _xy_ 2
#define _yy_ 3
#define _xz_ 4
#define _yz_ 5
#define _zz_ 6
#define _yx_ _xy_
#define _zx_ _xz_
#define _zy_ _yz_

#define _xxx_ 1
#define _xxy_ 2
#define _xxz_ 3
#define _xyy_ 4
#define _xyz_ 5
#define _xzz_ 6
#define _yyy_ 7
#define _yyz_ 8
#define _yzz_ 9
#define _zzz_ 10
#define _xyx_ _xxy_
#define _xzx_ _xxz_
#define _xzy_ _xyz_
#define _yxx_ _xxy_
#define _yxy_ _xyy_
#define _yxz_ _xyz_
#define _yyx_ _xyy_
#define _yzx_ _xyz_
#define _yzy_ _yyz_
#define _zxx_ _xxz_
#define _zxy_ _xyz_
#define _zxz_ _xzz_
#define _zyx_ _xyz_
#define _zyy_ _yyz_
#define _zyz_ _yzz_
#define _zzx_ _xzz_
#define _zzy_ _yzz_
#endif
