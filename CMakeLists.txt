cmake_minimum_required(VERSION 3.22)

project(Deterministic_FVM_Code LANGUAGES Fortran)

include(CTest)

# The explicit Fortran interfaces use the conventional LP64 LAPACK ABI.
# Leave BLA_VENDOR unset so CMake may select an optimized implementation;
# callers can still request one explicitly with -DBLA_VENDOR=<vendor>.
set(BLA_SIZEOF_INTEGER 4)
# CMake 4.2 probes NVPL without QUIET when vendor selection is open.  Probe
# it quietly first and disable only the redundant probe when NVPL is absent.
if((NOT DEFINED BLA_VENDOR OR BLA_VENDOR STREQUAL "All" OR
    BLA_VENDOR STREQUAL "NVPL") AND
   NOT DEFINED CMAKE_DISABLE_FIND_PACKAGE_nvpl)
    find_package(nvpl QUIET CONFIG)
    if(NOT nvpl_FOUND)
        set(CMAKE_DISABLE_FIND_PACKAGE_nvpl TRUE)
        set(_FVM_DISABLED_MISSING_NVPL TRUE)
    endif()
endif()
find_package(LAPACK REQUIRED)
if(_FVM_DISABLED_MISSING_NVPL)
    unset(CMAKE_DISABLE_FIND_PACKAGE_nvpl)
endif()
unset(_FVM_DISABLED_MISSING_NVPL)

set(FVM_CORE_SOURCES
    src/variables/MOD_kinds.f90
    src/variables/MOD_constants.f90
    src/variables/Variables.f90
    src/IO/MOD_IO_routines.f90
    src/MOD_Recurrent_Functions.f90
    src/solvers/MOD_BiCGSTAB_Solver.f90
    src/solvers/MOD_GaussElimination.f90
    src/solvers/MOD_LAPACK_Interfaces.f90
    src/solvers/MOD_Eigvals_Solver.f90
    src/accelerations/MOD_CMFD_Package.f90
    src/accelerations/MOD_pCMFD_Package.f90
    src/accelerations/MOD_lpCMFD_Package.f90
    src/accelerations/MOD_HGCMFD_Package.f90
    src/accelerations/MOD_pHGCMFD_Package.f90
    src/meshers/carthesian/MOD_Quad_Mesh.f90
    src/meshers/triangular/MOD_TriEq_Mesh.f90
    src/meshers/triangular/MOD_TriIso_Mesh.f90
    src/IO/MOD_Input_Reader.f90
)

add_library(fvm_core STATIC ${FVM_CORE_SOURCES})
target_link_libraries(fvm_core PRIVATE LAPACK::LAPACK)

set(FVM_CORE_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/fvm_core")
set_target_properties(fvm_core PROPERTIES
    Fortran_MODULE_DIRECTORY "${FVM_CORE_MODULE_DIRECTORY}"
)
target_include_directories(fvm_core PUBLIC
    "$<BUILD_INTERFACE:${FVM_CORE_MODULE_DIRECTORY}>"
)

add_executable(Deterministic_FVM_Code src/main/Header.f90)
target_link_libraries(Deterministic_FVM_Code PRIVATE fvm_core)
set_target_properties(Deterministic_FVM_Code PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/application"
)

add_library(fvm_debug_checks STATIC
    src/__debug_routines/MOD_DebugChecks.f90
)
target_link_libraries(fvm_debug_checks PRIVATE fvm_core)

set(FVM_DEBUG_MODULE_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/modules/fvm_debug_checks"
)
set_target_properties(fvm_debug_checks PROPERTIES
    Fortran_MODULE_DIRECTORY "${FVM_DEBUG_MODULE_DIRECTORY}"
)
target_include_directories(fvm_debug_checks PUBLIC
    "$<BUILD_INTERFACE:${FVM_DEBUG_MODULE_DIRECTORY}>"
)

function(fvm_set_fortran_options target_name)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(${target_name} PRIVATE -ffree-line-length-none)
        target_compile_options(${target_name} PRIVATE
            $<$<CONFIG:Debug>:-fcheck=all>
            $<$<CONFIG:Debug>:-fbacktrace>
        )
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
        if(WIN32)
            target_compile_options(${target_name} PRIVATE /extend-source)
            target_compile_options(${target_name} PRIVATE
                $<$<CONFIG:Debug>:/check:all>
                $<$<CONFIG:Debug>:/traceback>
            )
        else()
            target_compile_options(${target_name} PRIVATE -extend-source)
            target_compile_options(${target_name} PRIVATE
                $<$<CONFIG:Debug>:-check=all>
                $<$<CONFIG:Debug>:-traceback>
            )
        endif()
    endif()
endfunction()

fvm_set_fortran_options(fvm_core)
fvm_set_fortran_options(fvm_debug_checks)
fvm_set_fortran_options(Deterministic_FVM_Code)

# The legacy readers use input/... paths relative to the process working
# directory. Stage the current inputs in the same layout in the build tree.
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/input")
foreach(input_name
        Options_Input.inp
        Material_Input.inp
        Geometry_Input.inp
        CMFD_Grid.inp
        HGCMFD_Grid.inp)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/input/${input_name}"
        "${CMAKE_CURRENT_BINARY_DIR}/input/${input_name}"
        COPYONLY
    )
endforeach()
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
