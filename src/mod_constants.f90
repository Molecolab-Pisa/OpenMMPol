#include <openmmpol_const.h>

module mod_constants
    !! Contains several constants that are usefoul around the code. There are 
    !! physical constants (there should be no duplicate around!), default 
    !! force-field parameters that are used when no other input is specified, 
    !! and internal constants that are used as easy-to-remember names of code
    !! constants.
    use iso_c_binding
    
    implicit none

#ifdef USE_I8
    integer(kind=c_int64_t), parameter :: ip = c_int64_t
#else
    integer(kind=c_int32_t), parameter :: ip = c_int32_t
#endif
    !! Required precision for integer type
    integer(ip), parameter :: rp = c_double !! Required precision for real type
    integer(ip), parameter :: lp = c_bool

    ! Physical constants
    real(rp), parameter :: angstrom2au = OMMP_FORT_ANG2AU
    !! Conversion factor from \(\AA\) to A.U.
    real(rp), parameter :: kcalmol2au = OMMP_FORT_KCALMOL2AU
    !! Conversion factor from kcal mol\(^{-1}\) to A.U.
    real(rp), parameter :: au2kcalmol = OMMP_FORT_AU2KCALMOL
    !! Conversion factor from A.U. to kcal mol\(^{-1}\)
    real(rp), parameter :: pi = 4.0 * atan(1.0)
    !! Value of \(\pi\)
    real(rp), parameter :: deg2rad = pi / 180.0
    !! Conversion from arc degree to radians
    real(rp), parameter :: rad2deg = 180.0 / pi
    !! Conversion from radians to arc degree

    ! Force Field parameters
    real(rp), parameter :: mscale_wang_al(4) = [0.0_rp, 0.0_rp, 1.0_rp, 1.0_rp]
    real(rp), parameter :: pscale_wang_al(4) = [0.0_rp, 0.0_rp, 1.0_rp, 1.0_rp]
    real(rp), parameter :: dscale_wang_al(4) = [0.0_rp, 0.0_rp, 1.0_rp, 1.0_rp]
    real(rp), parameter :: uscale_wang_al(4) = [0.0_rp, 0.0_rp, 1.0_rp, 1.0_rp]
    
    real(rp), parameter :: mscale_wang_dl(4) = [1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp]
    real(rp), parameter :: pscale_wang_dl(4) = [1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp]
    real(rp), parameter :: dscale_wang_dl(4) = [1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp]
    real(rp), parameter :: uscale_wang_dl(4) = [1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp]
    
    real(rp), parameter :: mscale_amoeba(4) = [0.0_rp, 0.0_rp, 0.4_rp, 0.8_rp]
    real(rp), parameter :: pscale_amoeba(4) = [0.0_rp, 0.0_rp, 1.0_rp, 1.0_rp]
    real(rp), parameter :: dscale_amoeba(4) = [0.0_rp, 1.0_rp, 1.0_rp, 1.0_rp]
    real(rp), parameter :: uscale_amoeba(4) = [1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp]
    real(rp), parameter :: pscale_intra_amoeba(4) = &
                                              [0.0_rp, 0.0_rp, 0.5_rp, 1.0_rp]

            
    real(rp), parameter :: thole_scale_wang_al = 2.5874_rp
    !! Default Thole scaling factor used in Wang-AL force field
    real(rp), parameter :: thole_scale_wang_dl = 2.0580_rp
    !! Default Thole scaling factor used in Wang-DL force field
    
    ! Code constants
    real(rp), parameter :: eps_rp = epsilon(0.0_rp) * 100
    !! Threshold used to compare REALS for queality
    integer(ip), parameter :: ommp_str_char_max = OMMP_STR_CHAR_MAX

    ! Flag handler
    integer(ip), parameter :: ommp_ff_amoeba = OMMP_FF_AMOEBA
    !! Amoeba force field id
    integer(ip), parameter :: ommp_ff_wang_al = OMMP_FF_WANG_AL
    !! Wang AL force field id
    integer(ip), parameter :: ommp_ff_wang_dl = OMMP_FF_WANG_DL
    !! Wang DL force field id
    integer(ip), parameter :: ommp_ff_amber = OMMP_FF_AMBER
    !! Amber-like force field
    integer(ip), parameter :: ommp_ff_unknown = OMMP_FF_UNKNOWN
    !! Unknown force field

    integer(ip), parameter :: ommp_solver_cg = OMMP_SOLVER_CG 
    !! Conjugate gradients solver id
    integer(ip), parameter :: ommp_solver_diis = OMMP_SOLVER_DIIS
    !! DIIS solver id
    integer(ip), parameter :: ommp_solver_inversion = OMMP_SOLVER_INVERSION
    !! Matrix inversion solver id
    integer(ip), parameter :: ommp_solver_default = OMMP_SOLVER_DEFAULT
    !! Default value for solver
    integer(ip), parameter :: ommp_solver_none = OMMP_SOLVER_NONE
    !! Placeholder equivalent to not passing the argument, mainly for C interfaces

    integer(ip), parameter :: ommp_matv_incore = OMMP_MATV_INCORE
    !! Build matrix in memory to perform vector matrix multiplication
    !! in iterative solvers (uses dgemv with [[\mathcal O(2)]] memory usage
    integer(ip), parameter :: ommp_matv_direct = OMMP_MATV_DIRECT
    !! Vector matrix multiplication in iterative solvers are done in a direct
    !! fashion
    integer(ip), parameter :: ommp_matv_default = OMMP_MATV_DEFAULT
    !! Default value for matrix vector multiplication
    integer(ip), parameter :: ommp_matv_none = OMMP_MATV_NONE
    !! Placeholder equivalent to not passing the argument, mainly for C interfaces

    integer(ip), parameter :: ommp_verbose_debug = OMMP_VERBOSE_DEBUG
    !! Maximum verbosity level allowed
    integer(ip), parameter :: ommp_verbose_high = OMMP_VERBOSE_HIGH
    !! High level of verbosity
    integer(ip), parameter :: ommp_verbose_low = OMMP_VERBOSE_LOW
    !! Normal verbosity level
    integer(ip), parameter :: ommp_verbose_none = OMMP_VERBOSE_NONE
    !! All output should be suppressed at this level
    integer(ip), parameter :: ommp_verbose_default = OMMP_VERBOSE_DEFAULT
    !! All output should be suppressed at this level

    integer(ip), parameter :: AMOEBA_ROT_NONE = 0_ip 
    integer(ip), parameter :: AMOEBA_ROT_Z_THEN_X = 1_ip 
    integer(ip), parameter :: AMOEBA_ROT_BISECTOR = 2_ip
    integer(ip), parameter :: AMOEBA_ROT_Z_ONLY = 3_ip
    integer(ip), parameter :: AMOEBA_ROT_Z_BISECT = 4_ip
    integer(ip), parameter :: AMOEBA_ROT_3_FOLD = 5_ip 
    
    integer(ip), parameter :: OMMP_VDWTYPE_LJ = 0 
    integer(ip), parameter :: OMMP_VDWTYPE_BUF714 = 2
    integer(ip), parameter :: OMMP_RADRULE_ARITHMETIC = 1
    integer(ip), parameter :: OMMP_RADRULE_CUBIC = 2
    integer(ip), parameter :: OMMP_RADTYPE_RMIN = 1 
    integer(ip), parameter :: OMMP_EPSRULE_GEOMETRIC = 0 
    integer(ip), parameter :: OMMP_EPSRULE_HHG = 1

    integer(ip), parameter :: ommp_default_nl_sub = OMMP_DEFAULT_NL_SUB
    real(rp), parameter :: ommp_default_nl_cutoff = OMMP_DEFAULT_NL_CUTOFF
    
    integer(ip), parameter :: default_link_atom_n_eel_remove = OMMP_DEFAULT_LA_N_EEL_REMOVE
    real(rp), parameter :: default_link_atom_dist = OMMP_DEFAULT_LA_DIST

    integer(ip), parameter :: ommp_fmm_enable_thr = OMMP_FMM_ENABLE_THR
    integer(ip), parameter :: ommp_fmm_default_maxl_pol = OMMP_FMM_DEFAULT_MAXL_POL
    integer(ip), parameter :: ommp_fmm_default_maxl = OMMP_FMM_DEFAULT_MAXL
    real(rp), parameter :: ommp_fmm_default_cellsize = OMMP_FMM_DEFAULT_CELLSIZE
end module mod_constants
