module mod_constants
    !! Contains several constants that are usefoul around the code. There are 
    !! physical constants (there should be no duplicate around!), default 
    !! force-field parameters that are used when no other input is specified, 
    !! and internal constants that are used as easy-to-remember names of code
    !! constants.
    use mod_memory, only: ip, rp
    implicit none
  
    ! Physical constants
    real(rp), parameter :: angstrom2au = 1.8897261245650_rp
    !! Conversion factor from \(\AA\) to A.U.
    real(rp), parameter :: kcalmol2au = 1.59360109742136e-3_rp
    !! Conversion factor from kcal mol\(^{-1}\) to A.U.
    real(rp), parameter :: au2kcalmol = 627.5096080306_rp
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
    integer(ip), parameter :: OMMP_STR_CHAR_MAX = 1024

    ! Flag handler
    integer(ip), parameter :: OMMP_FF_AMOEBA = 1_ip
    !! Amoeba force field id
    integer(ip), parameter :: OMMP_FF_WANG_AL = 0_ip
    !! Wang AL force field id
    integer(ip), parameter :: OMMP_FF_WANG_DL = 0_ip
    !! Wang DL force field id

    integer(ip), parameter :: OMMP_SOLVER_CG = 1_ip 
    !! Conjugate gradients solver id
    integer(ip), parameter :: OMMP_SOLVER_DIIS = 2_ip 
    !! DIIS solver id
    integer(ip), parameter :: OMMP_SOLVER_INVERSION = 3_ip
    !! Matrix inversion solver id
    integer(ip), parameter :: OMMP_SOLVER_DEFAULT = OMMP_SOLVER_CG
    !! Default value for solver

    integer(ip), parameter :: OMMP_MATV_INCORE = 1_ip
    !! Build matrix in memory to perform vector matrix multiplication
    !! in iterative solvers (uses dgemv with [[\mathcal O(2)]] memory usage
    integer(ip), parameter :: OMMP_MATV_DIRECT = 2_ip
    !! Vector matrix multiplication in iterative solvers are done in a direct
    !! fashion
    integer(ip), parameter :: OMMP_MATV_DEFAULT = OMMP_MATV_DIRECT
    !! Default value for matrix vector multiplication

    integer(ip), parameter :: OMMP_AMOEBA_D = 1
    !! Index of direct (D) field and dipoles in AMOEBA FF (true dipoles)
    integer(ip), parameter :: OMMP_AMOEBA_P = 2
    !! Index of polarization (P) field and dipoles in AMOEBA FF (Lagrange multiplier)

    integer(ip), parameter :: OMMP_VERBOSE_DEBUG = 3_ip
    !! Maximum verbosity level allowed
    integer(ip), parameter :: OMMP_VERBOSE_HIGH = 2_ip
    !! High level of verbosity
    integer(ip), parameter :: OMMP_VERBOSE_LOW = 1_ip
    !! Normal verbosity level
    integer(ip), parameter :: OMMP_VERBOSE_NONE = 0_ip
    !! All output should be suppressed at this level

    integer(ip), parameter :: AMOEBA_ROT_NONE = 0_ip 
    integer(ip), parameter :: AMOEBA_ROT_Z_THEN_X = 1_ip 
    integer(ip), parameter :: AMOEBA_ROT_BISECTOR = 2_ip
    integer(ip), parameter :: AMOEBA_ROT_Z_ONLY = 3_ip
    integer(ip), parameter :: AMOEBA_ROT_Z_BISECT = 4_ip
    integer(ip), parameter :: AMOEBA_ROT_3_FOLD = 5_ip 

end module mod_constants
