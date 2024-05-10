set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs
  "${dir}/mod_adjacency_mat.f90"
  "${dir}/mod_bonded.f90"
  "${dir}/mod_c_interface.f90"
  "${dir}/mod_constants.F90"
  "${dir}/mod_electrostatics.F90"
  "${dir}/mod_geomgrad.F90"
  "${dir}/mod_inputloader.f90"
  "${dir}/mod_interface.F90"
  "${dir}/mod_io.F90"
  "${dir}/mod_iohdf5.F90"
  "${dir}/mod_jacobian_mat.f90"
  "${dir}/mod_link_atom.F90"
  "${dir}/mod_memory.f90"
  "${dir}/mod_mmpol.f90"
  "${dir}/mod_neighbors_list.F90"
  "${dir}/mod_nonbonded.f90"
  "${dir}/mod_polarization.F90"
  "${dir}/mod_prm.f90"
  "${dir}/mod_profiling.F90"
  "${dir}/mod_qm_helper.f90"
  "${dir}/mod_solvers.f90"
  "${dir}/mod_topology.f90"
  "${dir}/mod_utils.f90"
  "${dir}/ommp_pp.c"
  "${dir}/ommp_xyz2mmp.c"
  "${dir}/prm_keywords.f90"
  "${dir}/rotate_multipoles.F90"
  "${dir}/smartinput.c"
)

set_source_files_properties(
  "${dir}/mod_constants.F90"
  "${dir}/mod_electrostatics.F90"
  "${dir}/mod_geomgrad.F90"
  "${dir}/mod_interface.F90"
  "${dir}/mod_io.F90"
  "${dir}/mod_iohdf5.F90"
  "${dir}/mod_linkatom.F90"
  "${dir}/mod_neighbors_list.F90"
  "${dir}/mod_polarization.F90"
  "${dir}/mod_profiling.F90"
  "${dir}/rotate_multipoles.F90"
  PROPERTIES Fortran_PREPROCESS ON
)

set(srcs "${srcs}" PARENT_SCOPE)
