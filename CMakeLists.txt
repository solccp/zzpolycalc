if (NOT DEFINED CC)
  find_program (INTEL_C_COMPILER icc)
  if (INTEL_C_COMPILER)
    find_program(CMAKE_C_COMPILER NAMES $ENV{CC} icc PATHS ENV PATH NO_DEFAULT_PATH)
  endif()
  find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
endif()
if (NOT DEFINED FC)
  find_program (INTEL_FORTRAN_COMPILER ifort)
  if (INTEL_FORTRAN_COMPILER)
    find_program(CMAKE_Fortran_COMPILER NAMES $ENV{FC} ifort PATHS ENV PATH NO_DEFAULT_PATH)
  else()
    find_program(CMAKE_Fortran_COMPILER NAMES $ENV{FC} gfortran PATHS ENV PATH NO_DEFAULT_PATH)
  endif()
endif()


# Project name is not mandatory, but you should use it
project(ZZPolyCalc Fortran C)

cmake_minimum_required(VERSION 2.8.12..3.24.2)

# The macro below forces the build directory to be different from source directory:
#include(MacroOutOfSourceBuild)
 
#macro_ensure_out_of_source_build("${PROJECT_NAME} requires an out of source build.")

option(USE_SHA256 "Enable SHA256" OFF)





set(CMAKE_INSTALL_PREFIX .)

#set( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
#SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE) # do not strip non-standard libraries when make install

if (USE_SHA256 STREQUAL ON)
    find_library(SSL_LIBRARY NAMES crypto)
    if(SSL_LIBRARY)
      message("Found libcrypto. Will use this SSL library for SHA256")
    else()
      message("FATAL_ERROR libssl required for SHA256. Edit CMakeLists.txt and turn off SHA256 or specify SSL location via CMAKE_PREFIX_PATH")
    endif()
else()
   find_path(XXHASH_INCLUDE_DIR NAMES xxhash.h)
   if(NOT XXHASH_INCLUDE_DIR)
     message(WARNING "Didn't find xxhash.h. Will use slower MD5")
   endif()
   find_library(XXHASH_LIBRARY NAMES xxhash)
   if(NOT XXHASH_LIBRARY)
     message(WARNING "Didn't find XXHASH library. Will use slower MD5")
     find_library(SSL_LIBRARY NAMES crypto)
     if(SSL_LIBRARY)
       message("Found libcrypto. Will use this SSL library for MD5")
     elseif (USE_SHA256 STREQUAL ON)
       message("FATAL_ERROR libssl required for SHA256")
      endif()
  else()
     message("Will use XXHASH128")
     option(USE_XXHASH "enable XXHASH" ON)
  endif()
endif()


    


add_subdirectory(src)
enable_testing()
add_subdirectory(test)
