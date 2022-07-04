module mod_constants
    use mod_memory, only: ip, rp
    implicit none
  
    real(rp), parameter :: zero = 0.0_rp, pt5 = 0.5_rp, one = 1.0_rp, two = 2.0_rp, three = 3.0_rp, &
                              four = 4.0_rp, five = 5.0_rp, six = 6.0_rp, seven = 7.0_rp, nine = 9.0_rp, &
                              ten = 10.0_rp, f15 = 15.0_rp, f105 = 105.0_rp
    ! Physical constants
    real(rp), parameter :: angstrom2au = 1.8897261245650_rp 
    !! Conversion factor Angstrom -> Bohr

    ! Force Field parameters
    real(rp), parameter :: mscale_wang_al(4) = (/0.0_rp, 0.0_rp, 1.0_rp, 1.0_rp/)
    real(rp), parameter :: pscale_wang_al(5) = (/0.0_rp, 0.0_rp, 1.0_rp, 1.0_rp, 1.0_rp/)
    real(rp), parameter :: dscale_wang_al(4) = (/0.0_rp, 0.0_rp, 1.0_rp, 1.0_rp/)
    real(rp), parameter :: uscale_wang_al(4) = (/0.0_rp, 0.0_rp, 1.0_rp, 1.0_rp/)
    
    real(rp), parameter :: mscale_wang_dl(4) = (/1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp/)
    real(rp), parameter :: pscale_wang_dl(5) = (/1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp/)
    real(rp), parameter :: dscale_wang_dl(4) = (/1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp/)
    real(rp), parameter :: uscale_wang_dl(4) = (/1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp/)
    
    real(rp), parameter :: mscale_amoeba(4) = (/0.0_rp, 0.0_rp, 0.4_rp, 0.8_rp/)
    real(rp), parameter :: pscale_amoeba(5) = (/0.0_rp, 0.0_rp, 1.0_rp, 1.0_rp, 0.5_rp/)
    real(rp), parameter :: dscale_amoeba(4) = (/0.0_rp, 1.0_rp, 1.0_rp, 1.0_rp/)
    real(rp), parameter :: uscale_amoeba(4) = (/1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp/)
            
    real(rp),   parameter   :: a_wal = 2.5874_rp, a_wdl = 2.0580_rp
    
    ! Code constants
    real(rp), parameter :: eps_rp = epsilon(0.0_rp) * 100
    !! Threshold used to compare REALS for queality

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
    !! Index of direct (D) field and dipoles in AMOEBA FF
    integer(ip), parameter :: OMMP_AMOEBA_P = 2
    !! Index of polarization (P) field and dipoles in AMOEBA FF

    integer(ip), parameter :: OMMP_VERBOSE_DEBUG = 3_ip
    integer(ip), parameter :: OMMP_VERBOSE_HIGH = 2_ip
    integer(ip), parameter :: OMMP_VERBOSE_LOW = 1_ip
    integer(ip), parameter :: OMMP_VERBOSE_NONE = 0_ip

    integer(ip), parameter :: AMOEBA_ROT_NONE = 0_ip 
    integer(ip), parameter :: AMOEBA_ROT_Z_THEN_X = 1_ip 
    integer(ip), parameter :: AMOEBA_ROT_BISECTOR = 2_ip
    integer(ip), parameter :: AMOEBA_ROT_Z_ONLY = 3_ip
    integer(ip), parameter :: AMOEBA_ROT_Z_BISECT = 4_ip
    integer(ip), parameter :: AMOEBA_ROT_3_FOLD = 5_ip 

end module mod_constants
