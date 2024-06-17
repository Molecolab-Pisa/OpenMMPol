# Borrowed from https://github.com/dftbplus/dftbplus/blob/main/cmake/DftbPlusUtils.cmake

# Loads build settings
macro (get_config_arguments)

# Find custom config file; if not found, use global
  if(NOT DEFINED BUILD_CONFIG_FILE)
	  if(DEFINED ENV{OPENMMPOL_BUILD_CONFIG_FILE}
			  AND NOT "$ENV{OPENMMPOL_BUILD_CONFIG_FILE}" STREQUAL "")
		  set(BUILD_CONFIG_FILE "$ENV{OPENMMPOL_BUILD_CONFIG_FILE}")
    else()
      set(BUILD_CONFIG_FILE "${OMMP_ROOT_DIR}/config.cmake")
    endif()
  endif()
  message(STATUS "Reading global build config file: ${BUILD_CONFIG_FILE}")
  include(${BUILD_CONFIG_FILE})

endmacro()

# Sets up the build type.
function (setup_build_type)

  set(default_build_type "RelWithDebInfo")
  get_property(_multiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
  if(_multiConfig)
    set(CMAKE_CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo;Coverage")
    message(STATUS "Build type: Multi-Config (build type selected at the build step)")
  else()
    if(NOT CMAKE_BUILD_TYPE)
      message(STATUS "Build type: ${default_build_type} (default single-config)")
      set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Build type" FORCE)
      set_property(CACHE CMAKE_BUILD_TYPE PROPERTY HELPSTRING "Choose the type of build")
      set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo" "Coverage")
    else()
      message(STATUS "Build type: ${CMAKE_BUILD_TYPE} (manually selected single-config)")
    endif()
  endif()

endfunction()


# Tries to guess which toolchain to load based on the environment.
#
# Args:
#     toolchain [out]: Name of the selected toolchain or undefined if it could not be selected
#
function(guess_toolchain toolchain)

  if("${CMAKE_Fortran_COMPILER_ID}|${CMAKE_C_COMPILER_ID}" STREQUAL "GNU|GNU")
    set(_toolchain "gnu")
  elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel|IntelLLVM"
      AND "${CMAKE_C_COMPILER_ID}" MATCHES "Intel|IntelLLVM")
    set(_toolchain "intel")
  # elseif("${CMAKE_Fortran_COMPILER_ID}|${CMAKE_C_COMPILER_ID}" STREQUAL "NAG|GNU")
    # set(_toolchain "nag")
  # TODO: add toolchains for other compilers, e.g. intel
  else()
    set(_toolchain "generic")
  endif()

  set(${toolchain} "${_toolchain}" PARENT_SCOPE)

endfunction()


# Loads toolchain settings.
#
macro(load_toolchain_settings)

  if(NOT DEFINED TOOLCHAIN_FILE AND NOT "$ENV{OPENMMPOL_TOOLCHAIN_FILE}" STREQUAL "")
    set(TOOLCHAIN_FILE "$ENV{OPENMMPOL_TOOLCHAIN_FILE}")
  endif()
  if(NOT DEFINED TOOLCHAIN AND NOT "$ENV{OPENMMPOL_TOOLCHAIN}" STREQUAL "")
    set(TOOLCHAIN "$ENV{OPENMMPOL_TOOLCHAIN}")
  endif()
  if(NOT DEFINED TOOLCHAIN_FILE OR TOOLCHAIN_FILE STREQUAL "")
    if(NOT DEFINED TOOLCHAIN OR TOOLCHAIN STREQUAL "")
      guess_toolchain(TOOLCHAIN)
    endif()
    set(TOOLCHAIN_FILE ${CMAKE_CURRENT_SOURCE_DIR}/toolchain/${TOOLCHAIN}.cmake)
  endif()
  message(STATUS "Reading build environment specific toolchain file: ${TOOLCHAIN_FILE}")
  include(${TOOLCHAIN_FILE})
endmacro()


# Sets up the global compiler flags
#
macro(setup_global_compiler_flags)

  if(CMAKE_BUILD_TYPE)
    set(_buildtypes ${CMAKE_BUILD_TYPE})
  else()
    set(_buildtypes ${CMAKE_CONFIGURATION_TYPES})
  endif()
  foreach(_buildtype IN LISTS _buildtypes)
    foreach (lang IN ITEMS Fortran C CXX)
      string(TOUPPER "${_buildtype}" _buildtype_upper)
      set(CMAKE_${lang}_FLAGS " ${${lang}_FLAGS}")
      add_preprocessor_flags(CMAKE_${lang}_FLAGS)
      set(CMAKE_${lang}_FLAGS_${_buildtype_upper} " ${${lang}_FLAGS_${_buildtype_upper}}")

      message(STATUS "Flags for ${lang}-compiler (build type: ${_buildtype}): "
        "${CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${_buildtype_upper}}")
    endforeach()
  endforeach()
  unset(_buildtypes)
  unset(_buildtype)
  unset(_buildtype_upper)
endmacro()

function (add_preprocessor_flags flaglist)

  set(_flaglist "${${flaglist}}")

  if (WITH_HDF5)
    set(_flaglist  "-DUSE_HDF5 ${_flaglist}")
  endif()

  if (WITH_I8)
    set(_flaglist  "-DUSE_I8 ${_flaglist}")
  endif()

  set(${flaglist} ${_flaglist} PARENT_SCOPE)

endfunction()