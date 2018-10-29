add_executable(
    reduce_sum_0d
    ${CMAKE_CURRENT_SOURCE_DIR}/reduce/sum_0d.cpp)
add_executable(
    reduce_product_0d
    ${CMAKE_CURRENT_SOURCE_DIR}/reduce/product_0d.cpp)
add_executable(
    reduce_max_0d
    ${CMAKE_CURRENT_SOURCE_DIR}/reduce/max_0d.cpp)
add_executable(
    reduce_min_0d
    ${CMAKE_CURRENT_SOURCE_DIR}/reduce/min_0d.cpp)
add_executable(
    reduce_or_0d
    ${CMAKE_CURRENT_SOURCE_DIR}/reduce/or_0d.cpp)
add_executable(
    reduce_xor_0d
    ${CMAKE_CURRENT_SOURCE_DIR}/reduce/xor_0d.cpp)
add_executable(
    reduce_and_0d
    ${CMAKE_CURRENT_SOURCE_DIR}/reduce/and_0d.cpp)

add_executable(
    gather_0d
    ${CMAKE_CURRENT_SOURCE_DIR}/gather/gather_0d.cpp)
add_executable(
    gather_1d
    ${CMAKE_CURRENT_SOURCE_DIR}/gather/gather_1d.cpp)

add_executable(
    add_quant
    ${CMAKE_CURRENT_SOURCE_DIR}/register/add_quant.cpp)
add_executable(
    add_phase
    ${CMAKE_CURRENT_SOURCE_DIR}/register/add_phase.cpp)

set(TEST_TARGETS
    reduce_sum_0d;reduce_product_0d;reduce_max_0d;reduce_min_0d;reduce_or_0d;reduce_xor_0d;reduce_and_0d;gather_0d;gather_1d;add_quant;add_phase)

foreach (TEST_TARGET ${TEST_TARGETS})
    target_include_directories(${TEST_TARGET} PRIVATE ${PROJECT_SOURCE_DIR}/include)
    target_link_libraries(${TEST_TARGET} PRIVATE typhon)
    target_compile_features(${TEST_TARGET} PRIVATE cxx_std_11)

    add_test(
        NAME test_${TEST_TARGET}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ./${TEST_TARGET}
    )
endforeach (TEST_TARGET ${TEST_TARGETS})
