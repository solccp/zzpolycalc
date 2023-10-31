# Project name is not mandatory, but you should use it
project(ZZPolyCalc Fortran C)

cmake_minimum_required(VERSION 2.6)

# The macro below forces the build directory to be different from source directory:
#include(MacroOutOfSourceBuild)
 
#macro_ensure_out_of_source_build("${PROJECT_NAME} requires an out of source build.")

option(USE_SHA256 "enable SHA256" ON)

find_path(XXHASH_INCLUDE_DIR NAMES xxhash.h)
if(NOT XXHASH_INCLUDE_DIR)
  message(WARNING "Didn't find xxhash.h. Will use slower MD5")
endif()
  
find_library(XXHASH_LIBRARY NAMES xxhash)
if(NOT XXHASH_LIBRARY)
    message(WARNING "Didn't find XXHASH library. Will use slower MD5")
    find_library(SSL_LIBRARY NAMES crypto)
    if(SSL_LIBRARY)
      message("Found libcrypto. Will use this SSL library")
    endif()
else()
   option(USE_XXHASH "enable XXHASH" ON)
endif()


    
set(CMAKE_INSTALL_PREFIX .)

add_subdirectory(src)
