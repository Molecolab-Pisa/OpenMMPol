module mod_constants
    use mod_memory, only: ip, rp
    implicit none
  
    real(rp), parameter :: thres = 1.0e-8_rp
    real(rp), parameter :: zero = 0.0_rp, pt5 = 0.5_rp, one = 1.0_rp, two = 2.0_rp, three = 3.0_rp, &
                              four = 4.0_rp, five = 5.0_rp, six = 6.0_rp, seven = 7.0_rp, nine = 9.0_rp, &
                              ten = 10.0_rp, f15 = 15.0_rp, f105 = 105.0_rp

end module mod_constants
