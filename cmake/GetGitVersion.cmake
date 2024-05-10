
# Adapted from Chromium source code:
# https://chromium.googlesource.com/external/github.com/google/benchmark/+/refs/tags/v1.4.1/cmake/GetGitVersion.cmake
# 
#
# - Returns a version string from Git tags
#
# This function inspects the annotated git tags for the project and returns 
# CMake version format compliant string and an abbreviated commit hash as
# two variables.
#
#  get_git_version(<version_var> <commit_var>)
#
# - Example
#
# include(GetGitVersion)
# get_git_version(GIT_VERSION GIT_COMMIT)
#
# Requires CMake 2.8.11+

find_package(Git)
if(__get_git_version)
  return()
endif()
set(__get_git_version INCLUDED)
function(get_git_version var1 var2)
  if(GIT_EXECUTABLE)
      # git describe --tags | sed "s/\-/./g" | rev | cut -d"." -f2-  | rev | tr -d "\n"
      execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags
          COMMAND sed "s/-/./g"
          COMMAND rev
          COMMAND cut -d "." -f2- 
          COMMAND rev
          COMMAND tr -d "\n"
          # COMMAND sed "s/-/+r/;s/-/./g"
          WORKING_DIRECTORY ${OMMP_ROOT_DIR}
          RESULT_VARIABLE status
          OUTPUT_VARIABLE GIT_VERSION
          ERROR_QUIET)
      
      execute_process(COMMAND ${GIT_EXECUTABLE} log --pretty=format:'%h' -n 1
        COMMAND tr -d "\n"
        WORKING_DIRECTORY ${OMMP_ROOT_DIR}
        RESULT_VARIABLE status
        OUTPUT_VARIABLE GIT_COMMIT
        ERROR_QUIET
      )
        
  else()
      set(GIT_VERSION "0.0.0.0")
      set(GIT_COMMIT "notfound")
  endif()
  message("-- git Version: ${GIT_VERSION}, commit: ${GIT_COMMIT}")
  set(${var1} ${GIT_VERSION} PARENT_SCOPE)
  set(${var2} ${GIT_COMMIT} PARENT_SCOPE)
endfunction()