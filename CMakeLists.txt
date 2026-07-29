cmake_minimum_required(VERSION 3.16)

project(lbm_d2q9_cpp LANGUAGES CXX)

get_property(LBM_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT LBM_MULTI_CONFIG AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()
if(LBM_MULTI_CONFIG)
    set(CMAKE_CONFIGURATION_TYPES RelWithDebInfo Release Debug MinSizeRel CACHE STRING "" FORCE)
endif()

function(configure_lbm_output_directories)
    if(LBM_MULTI_CONFIG)
        foreach(config ${CMAKE_CONFIGURATION_TYPES})
            string(TOUPPER "${config}" config_upper)
            set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${config_upper}
                "${CMAKE_BINARY_DIR}/${config}/bin" PARENT_SCOPE)
            set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${config_upper}
                "${CMAKE_BINARY_DIR}/${config}/bin" PARENT_SCOPE)
            set(CMAKE_PDB_OUTPUT_DIRECTORY_${config_upper}
                "${CMAKE_BINARY_DIR}/${config}/bin" PARENT_SCOPE)
            set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${config_upper}
                "${CMAKE_BINARY_DIR}/${config}/lib" PARENT_SCOPE)
        endforeach()
    else()
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
            "${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/bin" PARENT_SCOPE)
        set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
            "${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/bin" PARENT_SCOPE)
        set(CMAKE_PDB_OUTPUT_DIRECTORY
            "${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/bin" PARENT_SCOPE)
        set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
            "${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/lib" PARENT_SCOPE)
    endif()
endfunction()

configure_lbm_output_directories()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(lbm
    src/collision.cpp
    src/collision3d.cpp
    src/geometry_mask.cpp
    src/geometry_mask3d.cpp
    src/solver.cpp
    src/solver3d.cpp
    src/two_phase_solver.cpp
    src/two_phase_solver3d.cpp
)

target_include_directories(lbm
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_definitions(lbm
    PUBLIC
        LBM_RESULT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/result"
        LBM_GEOMETRY_DIR="${CMAKE_CURRENT_SOURCE_DIR}/geometry"
)

function(configure_lbm_target target_name)
    if(MSVC)
        target_compile_options(${target_name} PRIVATE /FS /utf-8)
    endif()

    if(LBM_MULTI_CONFIG)
        foreach(config ${CMAKE_CONFIGURATION_TYPES})
            string(TOUPPER "${config}" config_upper)
            set_target_properties(${target_name}
                PROPERTIES
                    RUNTIME_OUTPUT_DIRECTORY_${config_upper} "${CMAKE_BINARY_DIR}/${config}/bin"
                    LIBRARY_OUTPUT_DIRECTORY_${config_upper} "${CMAKE_BINARY_DIR}/${config}/bin"
                    PDB_OUTPUT_DIRECTORY_${config_upper} "${CMAKE_BINARY_DIR}/${config}/bin"
                    COMPILE_PDB_OUTPUT_DIRECTORY_${config_upper} "${CMAKE_BINARY_DIR}/${config}/bin"
                    COMPILE_PDB_NAME_${config_upper} "${target_name}_compile"
                    ARCHIVE_OUTPUT_DIRECTORY_${config_upper} "${CMAKE_BINARY_DIR}/${config}/lib"
            )
        endforeach()
    else()
        set_target_properties(${target_name}
            PROPERTIES
                RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/bin"
                LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/bin"
                PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/bin"
                COMPILE_PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/bin"
                COMPILE_PDB_NAME "${target_name}_compile"
                ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/lib"
        )
    endif()

    if(MSVC)
        get_target_property(target_type ${target_name} TYPE)
        if(target_type STREQUAL "STATIC_LIBRARY")
            add_custom_command(TARGET ${target_name}
                POST_BUILD
                COMMAND "${CMAKE_COMMAND}"
                    -Dsrc="$<TARGET_FILE_DIR:${target_name}>/${target_name}.pdb"
                    -Ddst="${CMAKE_BINARY_DIR}/$<CONFIG>/bin/${target_name}.pdb"
                    -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/move_if_exists.cmake"
            )
        endif()
    endif()
endfunction()

if(MSVC)
    target_compile_options(lbm PRIVATE /W4 /permissive-)
else()
    target_compile_options(lbm PRIVATE -Wall -Wextra -Wpedantic)
endif()

configure_lbm_target(lbm)

add_executable(lbm_periodic examples/periodic_shear.cpp)
target_link_libraries(lbm_periodic PRIVATE lbm)
configure_lbm_target(lbm_periodic)

add_executable(lbm_cavity examples/lid_driven_cavity.cpp)
target_link_libraries(lbm_cavity PRIVATE lbm)
configure_lbm_target(lbm_cavity)

add_executable(lbm_capillary examples/capillary_displacement.cpp)
target_link_libraries(lbm_capillary PRIVATE lbm)
configure_lbm_target(lbm_capillary)

add_executable(lbm_droplet examples/droplet_impact.cpp)
target_link_libraries(lbm_droplet PRIVATE lbm)
configure_lbm_target(lbm_droplet)

add_executable(lbm_geometry_flow examples/geometry_flow.cpp)
target_link_libraries(lbm_geometry_flow PRIVATE lbm)
configure_lbm_target(lbm_geometry_flow)

add_executable(lbm_geometry_displacement examples/geometry_displacement.cpp)
target_link_libraries(lbm_geometry_displacement PRIVATE lbm)
configure_lbm_target(lbm_geometry_displacement)

add_executable(lbm_d3_cavity examples/d3_lid_driven_cavity.cpp)
target_link_libraries(lbm_d3_cavity PRIVATE lbm)
configure_lbm_target(lbm_d3_cavity)

add_executable(lbm_d3_geometry_flow examples/d3_geometry_flow.cpp)
target_link_libraries(lbm_d3_geometry_flow PRIVATE lbm)
configure_lbm_target(lbm_d3_geometry_flow)

add_executable(lbm_d3_geometry_displacement examples/d3_geometry_displacement.cpp)
target_link_libraries(lbm_d3_geometry_displacement PRIVATE lbm)
configure_lbm_target(lbm_d3_geometry_displacement)

add_executable(lbm_d3_droplet examples/d3_droplet_impact.cpp)
target_link_libraries(lbm_d3_droplet PRIVATE lbm)
configure_lbm_target(lbm_d3_droplet)

include(CTest)
if(BUILD_TESTING)
    add_executable(lbm_numerics_tests tests/lbm_numerics.cpp)
    target_link_libraries(lbm_numerics_tests PRIVATE lbm)
    configure_lbm_target(lbm_numerics_tests)
    add_test(NAME lbm_numerics COMMAND $<TARGET_FILE:lbm_numerics_tests>)

    add_executable(lbm_3d_numerics_tests tests/lbm_3d_numerics.cpp)
    target_link_libraries(lbm_3d_numerics_tests PRIVATE lbm)
    configure_lbm_target(lbm_3d_numerics_tests)
    add_test(NAME lbm_3d_numerics COMMAND $<TARGET_FILE:lbm_3d_numerics_tests>)
endif()
