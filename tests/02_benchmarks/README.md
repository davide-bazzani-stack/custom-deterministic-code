# Eigenvalue solver benchmark

Configure and build the recommended release benchmark from the repository
root:

```bash
cmake --preset test-release
cmake --build --preset eigen-benchmark-release --parallel 2
```

Run it with one BLAS/OpenMP thread so that the measurement compares the
solver backends rather than nested thread policies:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  ./tests/build/release/tests/fvm_eigen_benchmark
```

The benchmark reports both first-solve and steady-state timings. The first
dense solve includes lazy LU factorization; the first vectorized solve includes
its workspace allocation. The median columns are measured after these warm-up
calls, so they represent repeated solves with cached LU/workspaces. Dense and
vectorized execution order alternates during the seven timed repetitions.

The executable stops with an error if the backends do not perform equal
Arnoldi work, disagree numerically, or fail an independently recomputed true
residual. `dense_over_vector_speedup` is the steady-state dense median divided
by the vectorized median. Setup timings are single samples, exclude the first
solve, and are diagnostic only. For MKL, BLIS, or another threaded BLAS, also
set that vendor's thread-count variable to one.
