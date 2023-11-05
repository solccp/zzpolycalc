

if (TEST_XML)
    execute_process(COMMAND ${TEST_PROG} -Q ${SOURCEDIR}/${TEST_NAME}.xyz TIMEOUT 120 RESULT_VARIABLE HAD_ERROR OUTPUT_FILE output.txt ERROR_FILE error.txt) 
elseif (TEST_SMALLBUCK)
    execute_process(COMMAND ${TEST_PROG} -s 100 ${SOURCEDIR}/${TEST_NAME}.xyz TIMEOUT 120 RESULT_VARIABLE HAD_ERROR OUTPUT_FILE output.txt ERROR_FILE error.txt)
else ()
    execute_process(COMMAND ${TEST_PROG} ${SOURCEDIR}/${TEST_NAME}.xyz TIMEOUT 120 RESULT_VARIABLE HAD_ERROR OUTPUT_FILE output.txt ERROR_FILE error.txt) 
endif()

if(HAD_ERROR)
    message(FATAL_ERROR "Test failed")
endif()

if (TEST_XML)
    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files output.txt ${SOURCEDIR}/${TEST_NAME}.xml RESULT_VARIABLE DIFFERENT1)
else ()
    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files output.txt ${SOURCEDIR}/${TEST_NAME}.out RESULT_VARIABLE DIFFERENT1)
endif()
if(DIFFERENT1)
        message(FATAL_ERROR "Test failed - ZZ_polynomial files differ")
endif()



