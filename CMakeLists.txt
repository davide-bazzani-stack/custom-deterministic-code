cmake_minimum_required(VERSION 3.20)

project(Deterministic_FVM_Code LANGUAGES Fortran)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/output")

add_executable(Deterministic_FVM_Code
    MOD_IO_routines.f90
    Variables.f90
    MOD_Recurrent_Functions.f90
    MOD_BiCGSTAB_Solvers.f90
    MOD_GaussElimination.f90
    MOD_EigenvaluesSpectrumSolver.f90
    MOD_CMFD_Package.f90
    MOD_pCMFD_Package.f90
    MOD_lpCMFD_Package.f90
    MOD_HGCMFD_Package.f90
    MOD_pHGCMFD_Package.f90
    MOD_Quad_Mesh.f90
    MOD_TriEq_Mesh.f90
    MOD_TriIso_Mesh.f90
    MOD_DebugChecks.f90
    MOD_Input_Reader.f90
    Header.f90
)

set_target_properties(Deterministic_FVM_Code PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
)

# Several source lines exceed Fortran's default 132-column free-form limit.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(Deterministic_FVM_Code PRIVATE -ffree-line-length-none)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
    if(WIN32)
        target_compile_options(Deterministic_FVM_Code PRIVATE /extend-source)
    else()
        target_compile_options(Deterministic_FVM_Code PRIVATE -extend-source)
    endif()
endif()
