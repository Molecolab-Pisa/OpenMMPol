set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs
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
  "${dir}/ommp_pp.c"
  "${dir}/ommp_xyz2mmp.c"
  "${dir}/rotate_multipoles.F90"
  "${dir}/smartinput.c"
)

# list(
#   APPEND srcs_preprocess
#   "${dir}/mod_c_interface.F90"
#   "${dir}/mod_constants.F90"
#   "${dir}/mod_electrostatics.F90"
#   "${dir}/mod_geomgrad.F90"
#   "${dir}/mod_interface.F90"
#   "${dir}/mod_io.F90"
#   "${dir}/mod_iohdf5.F90"
#   "${dir}/mod_link_atom.F90"
#   "${dir}/mod_memory.F90"
#   "${dir}/mod_neighbors_list.F90"
#   "${dir}/mod_polarization.F90"
#   "${dir}/mod_prm.F90"
#   "${dir}/mod_profiling.F90"
#   "${dir}/rotate_multipoles.F90"
# )

set(srcs "${srcs}" PARENT_SCOPE)
# set(srcs_preprocess "${srcs_preprocess}" PARENT_SCOPE)

