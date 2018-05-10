add_executable(
    exchange_cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/exchange.cpp)

target_include_directories(exchange_cpp PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(exchange_cpp PRIVATE typhon)
