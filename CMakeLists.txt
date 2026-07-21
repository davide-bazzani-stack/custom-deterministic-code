cmake_minimum_required(VERSION 3.20)

project(Deterministic_FVM_Code LANGUAGES Fortran)

include(CTest)

set(FVM_CORE_SOURCES
    src/MOD_IO_routines.f90
    src/Variables.f90
    src/MOD_Recurrent_Functions.f90
    src/MOD_BiCGSTAB_Solvers.f90
    src/MOD_GaussElimination.f90
    src/MOD_EigenvaluesSpectrumSolver.f90
    src/MOD_CMFD_Package.f90
    src/MOD_pCMFD_Package.f90
    src/MOD_lpCMFD_Package.f90
    src/MOD_HGCMFD_Package.f90
    src/MOD_pHGCMFD_Package.f90
    src/MOD_Quad_Mesh.f90
    src/MOD_TriEq_Mesh.f90
    src/MOD_TriIso_Mesh.f90
    src/MOD_DebugChecks.f90
    src/MOD_Input_Reader.f90
)

add_library(fvm_core STATIC ${FVM_CORE_SOURCES})

set(FVM_CORE_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/fvm_core")
set_target_properties(fvm_core PROPERTIES
    Fortran_MODULE_DIRECTORY "${FVM_CORE_MODULE_DIRECTORY}"
)
target_include_directories(fvm_core PUBLIC
    "$<BUILD_INTERFACE:${FVM_CORE_MODULE_DIRECTORY}>"
)

add_executable(Deterministic_FVM_Code src/Header.f90)
target_link_libraries(Deterministic_FVM_Code PRIVATE fvm_core)
set_target_properties(Deterministic_FVM_Code PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/application"
)

function(fvm_set_fortran_options target_name)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(${target_name} PRIVATE -ffree-line-length-none)
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
        if(WIN32)
            target_compile_options(${target_name} PRIVATE /extend-source)
        else()
            target_compile_options(${target_name} PRIVATE -extend-source)
        endif()
    endif()
endfunction()

fvm_set_fortran_options(fvm_core)
fvm_set_fortran_options(Deterministic_FVM_Code)

# The legacy readers use fixed filenames relative to the process working
# directory. Stage the current inputs in the build tree without changing them.
foreach(input_name
        Options_Input.inp
        Material_Input.inp
        Geometry_Input.inp
        CMFD_Grid.inp
        HGCMFD_Grid.inp)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/input/${input_name}"
        "${CMAKE_CURRENT_BINARY_DIR}/${input_name}"
        COPYONLY
    )
endforeach()
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")

if(BUILD_TESTING AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt")
    add_subdirectory(tests)
endif()
