set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs
  "${dir}/fast_multipoles/mod_fmm.F90"
  "${dir}/fast_multipoles/mod_fmm_interface.F90"
  "${dir}/fast_multipoles/mod_fmm_utils.F90"
  "${dir}/fast_multipoles/mod_harmonics.F90"
  "${dir}/fast_multipoles/mod_octatree.F90"
  "${dir}/fast_multipoles/mod_ribtree.F90"
  "${dir}/fast_multipoles/mod_tree.F90"
  "${dir}/mod_iohdf5.F90"
  "${dir}/mod_adjacency_mat.F90"
  "${dir}/mod_bonded.F90"
  "${dir}/mod_c_interface.F90"
  "${dir}/mod_constants.F90"
  "${dir}/mod_electrostatics.F90"
  "${dir}/mod_geomgrad.F90"
  "${dir}/mod_inputloader.F90"
  "${dir}/mod_interface.F90"
  "${dir}/mod_io.F90"
  "${dir}/mod_jacobian_mat.F90"
  "${dir}/mod_link_atom.F90"
  "${dir}/mod_memory.F90"
  "${dir}/mod_mmpol.F90"
  "${dir}/mod_neighbors_list.F90"
  "${dir}/mod_nonbonded.F90"
  "${dir}/mod_polarization.F90"
  "${dir}/mod_prm.F90"
  "${dir}/mod_profiling.F90"
  "${dir}/mod_qm_helper.F90"
  "${dir}/mod_solvers.F90"
  "${dir}/mod_topology.F90"
  "${dir}/mod_utils.F90"
  "${dir}/rotate_multipoles.F90"
  "${dir}/smartinput.c"
)

set(srcs "${srcs}" PARENT_SCOPE)

