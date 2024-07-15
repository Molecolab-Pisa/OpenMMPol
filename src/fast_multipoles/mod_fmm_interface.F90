module fmmlib_interface
    use mod_fmm, only: fmm_type, &
                       fmm_tree_type, &
                       fmm_init, free_fmm, &
                       tree_p2m, tree_m2m, tree_m2l, tree_l2l, &
                       fmm_solve, &
                       cart_prop_at_ipart, cart_propfar_at_ipart, cart_propnear_at_ipart
    use mod_tree, only: free_tree
    use mod_ribtree, only: init_as_ribtree
    use mod_octatree, only: init_as_octatree
    use mod_constants, only: fmmlib_real => rp, &
                             fmmlib_int => ip, &
                             fmmlib_pi => pi
    use mod_harmonics, only: fmm_m2l
    use mod_fmm_utils, only: fmm_error

    implicit none

end module
