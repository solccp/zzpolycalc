file(GLOB TEST_FILES *.xyz)
foreach(loop_var IN LISTS TEST_FILES)
    get_filename_component(FNAME ${loop_var} NAME_WE)
    add_test(NAME ${FNAME}
        COMMAND ${CMAKE_COMMAND}
        -DTEST_PROG=$<TARGET_FILE:ZZPolyCalc>
        -DTEST_NAME=${FNAME}
        -DSOURCEDIR=${CMAKE_HOME_DIRECTORY}/test
        -P ${CMAKE_HOME_DIRECTORY}/cmake/script/runtest.cmake)
endforeach()
foreach(loop_var IN LISTS TEST_FILES)
    get_filename_component(FNAME ${loop_var} NAME_WE)
    add_test(NAME ${FNAME}_smallbuck
        COMMAND ${CMAKE_COMMAND}
        -DTEST_PROG=$<TARGET_FILE:ZZPolyCalc>
        -DTEST_NAME=${FNAME}
        -DTEST_SMALLBUCK=ON
        -DSOURCEDIR=${CMAKE_HOME_DIRECTORY}/test
        -P ${CMAKE_HOME_DIRECTORY}/cmake/script/runtest.cmake)
endforeach()
file(GLOB TEST_FILES *.con)
foreach(loop_var IN LISTS TEST_FILES)
    get_filename_component(FNAME ${loop_var} NAME_WE)
    add_test(NAME ${FNAME}
        COMMAND ${CMAKE_COMMAND}
        -DTEST_PROG=$<TARGET_FILE:ZZPolyCalc>
        -DTEST_NAME=${FNAME}
        -DTEST_ADJ=ON        
        -DSOURCEDIR=${CMAKE_HOME_DIRECTORY}/test
        -P ${CMAKE_HOME_DIRECTORY}/cmake/script/runtest.cmake)
endforeach()
