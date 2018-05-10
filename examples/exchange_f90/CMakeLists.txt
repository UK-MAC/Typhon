add_executable(
    exchange_f90
    ${CMAKE_CURRENT_SOURCE_DIR}/exchange.f90)

target_include_directories(exchange_f90 PRIVATE ${PROJECT_SOURCE_DIR}/build)
#target_link_libraries(exchange_f90 PRIVATE typhon_f90)
target_link_libraries(exchange_f90 PRIVATE typhon)
