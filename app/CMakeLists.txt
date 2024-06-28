add_executable(ommp_pp "${CMAKE_CURRENT_SOURCE_DIR}/ommp_pp.c")
target_link_libraries(ommp_pp openmmpol)

add_executable(ommp_xyz2mmp "${CMAKE_CURRENT_SOURCE_DIR}/ommp_xyz2mmp.c")
target_link_libraries(ommp_xyz2mmp openmmpol)
