
set(CMAKE_INCLUDE_CURRENT_DIR ON)
file(GLOB ZZPolyCalc_SOURCES *.F90 *.f90 *.f95)
if(XXHASH_INCLUDE_DIR)
  file(GLOB XXHASHWRAPPER xxhashwrapper.c)
  include_directories(${XXHASH_INCLUDE_DIR})
else()
  if(NOT SSL_LIBRARY)
    file(GLOB MD5 md5.c)
  endif()
endif()

if(USE_XXHASH STREQUAL ON)
    add_definitions( -DUSE_XXHASH )
endif()




add_executable(ZZPolyCalc ${ZZPolyCalc_SOURCES} ${MD5} ${XXHASHWRAPPER} )

#if( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )
#        set_target_properties(ZZPolyCalc PROPERTIES LINK_FLAGS " ")
#        set_target_properties(ZZPolyCalc PROPERTIES COMPILE_FLAGS "")
#endif()


if(XXHASH_LIBRARY)
  target_link_libraries(ZZPolyCalc ${XXHASH_LIBRARY})
endif()
if(SSL_LIBRARY)
  if(USE_SHA256 STREQUAL ON)
    add_definitions( -DUSE_SHA256 )
  endif()
  target_link_libraries(ZZPolyCalc ${SSL_LIBRARY})
endif()



INSTALL(TARGETS ZZPolyCalc
  RUNTIME DESTINATION ../bin
)
