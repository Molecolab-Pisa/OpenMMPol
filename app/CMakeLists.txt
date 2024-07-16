add_executable(ommp_xyz2mmp "${CMAKE_CURRENT_SOURCE_DIR}/ommp_xyz2mmp.c")
target_link_libraries(ommp_xyz2mmp openmmpol)

# CI custom targets; do not touch

if(WITH_HDF5)
    add_executable(ommp_pp "${CMAKE_CURRENT_SOURCE_DIR}/ommp_pp.c")
    target_link_libraries(ommp_pp openmmpol)
    add_custom_target(openmmpol_utils DEPENDS ommp_xyz2mmp ommp_pp)
else()
    add_custom_target(openmmpol_utils DEPENDS ommp_xyz2mmp)
endif()

