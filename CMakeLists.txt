cmake_minimum_required(VERSION 3.10)

if (NOT DEFINED ENABLE_FORTRAN)
    set(ENABLE_FORTRAN ON)
endif (NOT DEFINED ENABLE_FORTRAN)

if (ENABLE_FORTRAN)
    project(
        typhon
        VERSION 3.0.0
        DESCRIPTION "Typhon is a distributed communications library for unstructured mesh applications."
        LANGUAGES CXX Fortran
    )
else (ENABLE_FORTRAN)
    project(
        typhon
        VERSION 3.0.0
        DESCRIPTION "Typhon is a distributed communications library for unstructured mesh applications."
        LANGUAGES CXX
    )
endif (ENABLE_FORTRAN)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if (NOT DEFINED BUILD_TYPHON_DUMP)
    set(BUILD_TYPHON_DUMP OFF)
endif (NOT DEFINED BUILD_TYPHON_DUMP)

if (NOT DEFINED BUILD_DOCS)
    set(BUILD_DOCS OFF)
endif (NOT DEFINED BUILD_DOCS)


# ------------------------------------------------------------------------------
# Get VCS info
# ------------------------------------------------------------------------------
# Get the current working branch
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Add git info
add_definitions("-DTYPHON_GIT_BRANCH=${GIT_BRANCH}")
add_definitions("-DTYPHON_GIT_COMMIT_HASH=${GIT_COMMIT_HASH}")



# ------------------------------------------------------------------------------
# Shared library
# ------------------------------------------------------------------------------
set(TYPHON_SRC
    ${PROJECT_SOURCE_DIR}/src/collect.cpp
    ${PROJECT_SOURCE_DIR}/src/decomposition.cpp
    ${PROJECT_SOURCE_DIR}/src/dt_reduce.cpp
    ${PROJECT_SOURCE_DIR}/src/exchange.cpp
    ${PROJECT_SOURCE_DIR}/src/keys.cpp
    ${PROJECT_SOURCE_DIR}/src/quantity.cpp
    ${PROJECT_SOURCE_DIR}/src/register.cpp
    ${PROJECT_SOURCE_DIR}/src/schedule.cpp
    ${PROJECT_SOURCE_DIR}/src/core.cpp
    ${PROJECT_SOURCE_DIR}/src/utilities.cpp
    ${PROJECT_SOURCE_DIR}/src/serialise.cpp
    ${PROJECT_SOURCE_DIR}/src/distribute.cpp

    ${PROJECT_SOURCE_DIR}/src/distribute/ns.cpp
    ${PROJECT_SOURCE_DIR}/src/distribute/element_cloud.cpp
    ${PROJECT_SOURCE_DIR}/src/distribute/node_cloud.cpp
    ${PROJECT_SOURCE_DIR}/src/distribute/layer_info.cpp
)

add_library(typhon SHARED ${TYPHON_SRC})
set_target_properties(typhon PROPERTIES VERSION ${PROJECT_VERSION})

set (TYPHON_PUBLIC_HEADERS "include/typhon.h")
set_target_properties(typhon PROPERTIES PUBLIC_HEADER "${TYPHON_PUBLIC_HEADERS}")

target_compile_features(typhon PRIVATE cxx_std_11)
target_compile_options(typhon PRIVATE -Wall -Wextra -pedantic)

target_include_directories(typhon PRIVATE include)

if (CMAKE_BUILD_TYPE MATCHES "Debug")
    add_definitions("-DTYPHON_DEBUG")
endif (CMAKE_BUILD_TYPE MATCHES "Debug")



# ------------------------------------------------------------------------------
# Fortran module
# ------------------------------------------------------------------------------
if (ENABLE_FORTRAN)
    target_sources(typhon PRIVATE ${PROJECT_SOURCE_DIR}/src/typhon.f90)

    list (APPEND TYPHON_PUBLIC_HEADERS "build/typhon.mod")
    set_target_properties(typhon PROPERTIES PUBLIC_HEADER
        "${TYPHON_PUBLIC_HEADERS}")

    target_compile_options(typhon PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-std=f2008>)
    target_compile_options(typhon PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-c-binding-type>)
endif (ENABLE_FORTRAN)



# ------------------------------------------------------------------------------
# Utility binary
# ------------------------------------------------------------------------------
if (BUILD_TYPHON_DUMP)
    add_executable(
        typhon_dump
        ${PROJECT_SOURCE_DIR}/src/typhon_dump.cpp)

    set_target_properties(typhon_dump PROPERTIES VERSION ${PROJECT_VERSION})
    target_link_libraries(typhon_dump PRIVATE typhon)
    target_include_directories(typhon_dump PRIVATE include)
endif (BUILD_TYPHON_DUMP)



# ------------------------------------------------------------------------------
# MPI
# ------------------------------------------------------------------------------
set(MPI_CXX_SKIP_MPICXX "true")
if (ENABLE_FORTRAN)
    find_package(MPI COMPONENTS CXX Fortran REQUIRED)
else (ENABLE_FORTRAN)
    find_package(MPI COMPONENTS CXX REQUIRED)
endif (ENABLE_FORTRAN)

target_include_directories(typhon PUBLIC ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(typhon PUBLIC ${MPI_CXX_LIBRARIES})

if (ENABLE_FORTRAN)
    target_include_directories(typhon PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
    target_link_libraries(typhon PUBLIC ${MPI_Fortran_LIBRARIES})
endif (ENABLE_FORTRAN)

target_compile_options(typhon PUBLIC
    $<$<COMPILE_LANGUAGE:CXX>:${MPI_CXX_COMPILE_OPTIONS}>)
set_target_properties(typhon PROPERTIES LINK_FLAGS ${MPI_CXX_LINK_FLAGS})

if (ENABLE_FORTRAN)
    target_compile_options(typhon PUBLIC
        $<$<COMPILE_LANGUAGE:Fortran>:${MPI_Fortran_COMPILE_OPTIONS}>)
    set_target_properties(typhon PROPERTIES LINK_FLAGS
        ${MPI_Fortran_LINK_FLAGS})
endif (ENABLE_FORTRAN)



# ------------------------------------------------------------------------------
# yaml-cpp
# ------------------------------------------------------------------------------
if (BUILD_TYPHON_DUMP)
    find_package(YamlCpp)
    if (YamlCpp_FOUND)
        add_definitions("-DTYPHON_YAMLCPP_SUPPORT")
        target_include_directories(typhon_dump PRIVATE ${YamlCpp_INCLUDES})
        target_link_libraries(typhon_dump PRIVATE ${YamlCpp_LIBRARIES})

        target_sources(typhon_dump PRIVATE
            ${PROJECT_SOURCE_DIR}/src/yaml_output.cpp)
    endif (YamlCpp_FOUND)
endif (BUILD_TYPHON_DUMP)



# ------------------------------------------------------------------------------
# Doxygen
# ------------------------------------------------------------------------------
if (BUILD_DOCS)
    find_package(Doxygen REQUIRED)

    # set input and output files
    set(Doxygen_IN ${PROJECT_SOURCE_DIR}/doc/Doxyfile.in)
    set(Doxygen_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # request to configure the file
    configure_file(${Doxygen_IN} ${Doxygen_OUT} @ONLY)

    # Note the option ALL which allows to build the docs together with the application
    add_custom_target(doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${Doxygen_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)

endif (BUILD_DOCS)



# ------------------------------------------------------------------------------
# Installation
# ------------------------------------------------------------------------------
include(GNUInstallDirs)
install(
    TARGETS typhon
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

if (BUILD_TYPHON_DUMP)
    install(
        TARGETS typhon_dump
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif (BUILD_TYPHON_DUMP)


# ------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------
enable_testing()
add_subdirectory(test)



# ------------------------------------------------------------------------------
# Examples
# ------------------------------------------------------------------------------
add_subdirectory(examples)
