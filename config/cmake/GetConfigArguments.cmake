# Loads global build arguments
macro (GetConfigArguments)

# Find custom config file; if not found, use global
  if(NOT DEFINED BUILD_CONFIG_FILE)
	  if(DEFINED ENV{OPENMMPOL_BUILD_CONFIG_FILE}
			  AND NOT "$ENV{OPENMMPOL_BUILD_CONFIG_FILE}" STREQUAL "")
		  set(BUILD_CONFIG_FILE "$ENV{OPENMMPOL_BUILD_CONFIG_FILE}")
    else()
      set(BUILD_CONFIG_FILE "${CMAKE_CURRENT_SOURCE_DIR}/config.cmake")
    endif()
  endif()
  message(STATUS "Reading global build config file: ${BUILD_CONFIG_FILE}")
  include(${BUILD_CONFIG_FILE})

endmacro()