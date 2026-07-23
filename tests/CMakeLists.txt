set(TESTDRIVE_SOURCE_DIRECTORY
    "${CMAKE_CURRENT_SOURCE_DIR}/test-drive-main/test-drive-main/src"
)

add_executable(fvm_unit_tests
    "${TESTDRIVE_SOURCE_DIRECTORY}/testdrive_version.f90"
    "${TESTDRIVE_SOURCE_DIRECTORY}/testdrive.F90"
    00_unit/TEST_BiCGSTAB_hybrid_suite.f90
    00_unit/TEST_GaussElimination_suite.f90
    00_unit/TEST_Eigvals_Solver_suite.f90
    00_unit/MAIN_Test-Wrapper.f90
)

target_compile_definitions(fvm_unit_tests PRIVATE
    WITH_QP=0
    WITH_XDP=0
)
target_link_libraries(fvm_unit_tests PRIVATE fvm_core)
fvm_set_fortran_options(fvm_unit_tests)

set_target_properties(fvm_unit_tests PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/fvm_unit_tests"
)

add_test(NAME unit_solvers COMMAND fvm_unit_tests)
set_tests_properties(unit_solvers PROPERTIES
    LABELS "unit;solver"
    TIMEOUT 30
)

add_executable(fvm_eigen_benchmark
    02_benchmarks/BENCH_Eigvals_Solver.f90
)
target_link_libraries(fvm_eigen_benchmark PRIVATE fvm_core)
fvm_set_fortran_options(fvm_eigen_benchmark)
set_target_properties(fvm_eigen_benchmark PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/fvm_eigen_benchmark"
)
