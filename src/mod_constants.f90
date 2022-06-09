module mod_constants
    use mod_memory, only: ip, rp
    implicit none
  
    real(rp), parameter :: zero = 0.0_rp, pt5 = 0.5_rp, one = 1.0_rp, two = 2.0_rp, three = 3.0_rp, &
                              four = 4.0_rp, five = 5.0_rp, six = 6.0_rp, seven = 7.0_rp, nine = 9.0_rp, &
                              ten = 10.0_rp, f15 = 15.0_rp, f105 = 105.0_rp
    real(rp), parameter :: angstrom2au = 1.8897261245650_rp 
    !! Conversion factor Angstrom -> Bohr


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
    
    real(rp), parameter :: eps_rp = epsilon(0.0_rp) * 100
    real(rp), parameter :: thres = 1e-8
end module mod_constants
