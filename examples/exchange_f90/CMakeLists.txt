add_executable(
    exchange_f90
    ${CMAKE_CURRENT_SOURCE_DIR}/exchange.f90)

target_include_directories(exchange_f90 PRIVATE ${CMAKE_BINARY_DIR})
target_link_libraries(exchange_f90 PRIVATE typhon)
