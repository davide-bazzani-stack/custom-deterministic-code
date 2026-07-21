# Refactoring roadmap

This is the ordered work queue for the first refactoring. Detailed architecture,
defect evidence, tests, and invariants are recorded in `temp_refactoring.md`.

Priority meanings:

- **P0 — blocking:** establish tests, verify/fix solver correctness, and measure the
  highest-cost numerical paths before structural refactoring.
- **P1 — first refactoring:** establish controlled baselines and introduce the main
  ownership/module/strategy boundaries without changing input behavior.
- **P2 — follow-up quality:** complete lower-risk decomposition, numerical
  qualification, tooling, and maintainability work.
- **P3 — deferred input work:** explicitly excluded from the first refactoring.
- **P4 — cleanup/decisions:** repository hygiene and retirement decisions.

Rules for every implementation step:

- [ ] Do not mix a numerical behavior fix with a structural move.
- [ ] Run strict Debug tests before and after each step.
- [ ] Require residual/conservation/status checks, not compilation alone.
- [ ] Run cases that write files only in isolated copied working directories.
- [ ] Preserve input syntax and reader behavior throughout the first refactoring.
- [ ] Keep pure stateless kernels procedural; use OOP only for ownership, lifecycle,
  validation, and interchangeable strategies.

## P0 — Tests, linear-solver correctness, and measured optimization

### P0.1 Establish the safety net

- [ ] Establish an immutable source baseline before implementation.
  - [ ] Restore or initialize usable version-control metadata; the present `.git`
    directory is empty.
  - [ ] Record compiler/version, flags, and hashes of all input decks used for
    numerical baselines.
- [ ] Split the CMake build into a production `fvm_core` library and a thin
  executable without changing runtime behavior.
  - [ ] Link tests and the executable to the same library target.
  - [ ] Use explicit source lists; do not recursively glob the dormant `HGCMFD/`
    copy.
- [ ] Create the separate singular `test/` directory using the reference TestDrive
  collector/driver logic.
  - [ ] Vendor or pin one known `testdrive.F90` version and record its provenance.
  - [ ] Compile it with the required `WITH_QP=0` and `WITH_XDP=0` definitions.
  - [ ] Use one suite module per concern and one explicit wrapper that registers
    `new_testsuite` collectors, calls `run_testsuite`, accumulates failures, and
    exits nonzero.
  - [ ] Integrate CTest with labels/timeouts and `--output-on-failure` support.
  - [ ] Isolate Fortran `.mod` output by target/configuration.
  - [ ] Copy/configure fixtures into the test build tree; do not hard-code the
    reference repository path.
  - [ ] Disable parallel TestDrive execution for suites touching legacy fixed I/O
    units.
  - [ ] Add GNU and Intel strict-Debug flag profiles: warnings, bounds checks,
    traceback/backtrace, uninitialized diagnostics where available, and
    floating-point exception checks.
- [ ] Add benchmark targets separately from pass/fail CTest tests.

### P0.2 Characterize every linear operator and solver

- [ ] Add `TEST_stencil_matvec_suite`.
  - [ ] Compare Classic/Quad packed matvec with an independently assembled dense
    operator for arbitrary vectors and every basis vector.
  - [ ] Cover corners, edges, first/last row, rectangular grids, and smallest valid
    dimensions under bounds checking.
  - [ ] Compare TriIso packed matvec for all four triangle orientations with an
    independent neighbor-map operator.
  - [ ] Characterize EqTri offsets even though the fine-mesh implementation file is
    currently empty.
  - [ ] Compare HG arrow matvec with `matmul(to_dense(A), x)`.
- [ ] Add `TEST_bicgstab_suite` for every structured variant.
  - [ ] Exact initial solution.
  - [ ] Zero RHS/zero solution.
  - [ ] One-step diagonal solution.
  - [ ] Known nonsymmetric diagonally dominant systems and nonzero initial guesses.
  - [ ] Large/small scaling ranges.
  - [ ] `max_iter=0` and maximum-iteration failure.
  - [ ] Singular/zero operator and explicit rho/alpha/omega/t-norm breakdown cases.
  - [ ] Nonfinite-input handling.
  - [ ] On success require finite values and
    `||b-Ax||_2 <= atol + rtol*||b||_2`.
  - [ ] On failure require deterministic status/reason and no NaN, Inf, or hang.
- [ ] Add `TEST_dense_solver_suite` for `GE_Expl_Main`.
  - [ ] 1x1, identity, pivot-required, nonsymmetric well-conditioned, and
    fixed-seed diagonally dominant systems.
  - [ ] Exact singular and near-singular systems.
  - [ ] Scale-range cases.
  - [ ] Verify residual and known solution.
  - [ ] Verify input matrix/RHS remain unchanged when that is the API contract.
  - [ ] Require a defined singular/near-singular status.
- [ ] Add `TEST_arrow_solver_suite`.
  - [ ] Compare dense and vector/arrow migration assembly.
  - [ ] Compare dense reference, Schur complement, and arrow BiCGSTAB residuals.
  - [ ] Cover n=1, n=2, nominal, and production-like sizes.
  - [ ] Cover zero arrow diagonals and near-zero Schur complement.
- [ ] Add `TEST_solver_parity_suite` for equivalent dense, banded, and arrow
  systems.
- [ ] Reuse the reference testbed's small independently specified fixtures only
  after copying, sanitizing, and documenting them.

### P0.3 Repair confirmed correctness defects behind failing tests

- [ ] Correct the three banded BiCGSTAB implementations.
  - [ ] Test the initial residual before iteration.
  - [ ] Restore the valid `s`-residual early exit.
  - [ ] Add scale-aware breakdown/nonfinite checks for rho, alpha denominator,
    omega, and `dot(t,t)`.
  - [ ] Use residual-based success criteria.
  - [ ] Return status, reason, iterations, and final residual.
- [ ] Increment `iter_inn_lo` in both Quad and TriIso and report limit exhaustion
  as failure/non-convergence rather than success.
- [ ] Propagate linear, inner, and outer iteration failures to one top-level result.
- [ ] Ensure the returned extended flux corresponds to the final reduced flux,
  including acceleration on the last permitted outer iteration.
- [ ] Correct dense Gaussian elimination.
  - [ ] Replace arithmetic row swapping with a safe temporary/pivot representation.
  - [ ] Detect singular/near-singular pivots and return a status.
  - [ ] Use the authoritative real kind consistently.
- [ ] Replace every `product(Phi)>0` gate with a finite elementwise condition that
  expresses the intended physics.
- [ ] Correct `Accel_Ini` use of incoming `Opt_gl`; it currently declares
  `intent(out)` and then reads values populated by `Options_Card`.
- [ ] Make the missing `HGCMFD_Grid.inp` path deterministic; it currently leaves
  `box_mat_ID`/related mapping state undefined.
- [ ] Handle no-acceleration calls without passing an unallocated `Elem_CMFD` as an
  assumed-shape actual argument.
- [ ] Reject TriEq and lpCMFD explicitly until implementations exist.
- [ ] Add defined guards/status for zero normalization, homogenization, modulation,
  face-area/distance, correction-coefficient, and current-map denominators.
- [ ] Validate every `findloc` result before it is used as a subscript.

### P0.4 Benchmark and optimize the solvers

- [ ] Build reproducible benchmarks for small, medium, and production-like coarse
  systems; record residual parity, iterations, factorization count, median time,
  allocations, peak memory, compiler, and flags.
- [ ] Benchmark current repeated dense GE against factor-once/solve-many LU.
- [ ] Evaluate LAPACK `GETRF/GETRS`, `GESV`, and/or banded routines if the dependency
  is acceptable.
- [ ] Validate and benchmark the O(n) HG/pHG Schur-complement solver against dense
  GE and arrow BiCGSTAB.
- [ ] Reuse BiCGSTAB workspaces instead of allocating seven vectors per solve.
- [ ] Stop building both dense and structured matrices unless an explicit backend
  selection requires both.
- [ ] Store/factor group blocks rather than a full `(n_coarse*n_group)^2` dense
  matrix.
- [ ] Benchmark direct procedure calls against type-bound/abstract matvec dispatch
  before putting dynamic dispatch inside hot loops.
- [ ] Select production solver backends from correctness and benchmark evidence;
  keep dense solve as a reference/debug backend where appropriate.

## P1 — Controlled baselines and first modular/OOP refactoring

### P1.1 Establish trusted behavior baselines

- [ ] Add isolated end-to-end cases for Quad with no acceleration, CMFD, pCMFD,
  HGCMFD, and pHGCMFD where supported.
- [ ] Add corresponding TriIso cases.
- [ ] Run every supplied case directory in a copied temporary working directory.
- [ ] Add missing-grid, unsupported-mode, malformed-mapping, zero-power, and
  solver-failure integration cases.
- [ ] Record final eigenvalue, iteration counts, normalized residuals, flux/power
  norms, conservation metrics, and numerical tolerances.
- [ ] Do not bless the existing root or `output/` `.out` files without provenance;
  they currently contain conflicting result sets.

### P1.2 Introduce aggregate ownership without changing input behavior

- [ ] Add one `simulation_case_t` that owns options, materials/XS, fine mesh,
  optional coarse mesh, workspaces, and results.
- [ ] Add a `simulation_result_t` containing flux, eigenvalue, power, convergence
  history, timings, and final status.
- [ ] Keep `Options_Card`, `Material_Card`, `Geometry_Card`, and the coarse-card
  parser behavior behind a legacy compatibility adapter.
- [ ] Preserve existing flattened layouts during the first mechanical transition.
- [ ] Centralize full/reduced/group index calculations and validate array shapes.
- [ ] Do not create polymorphic per-cell objects; move toward contiguous arrays
  inside aggregate mesh/model types.

### P1.3 Establish focused module boundaries

- [ ] Create one authoritative `kinds_m` using `real64`; remove mixed `dp_var`,
  `doupr`, `real(8)`, and unintended default-real use incrementally.
- [ ] Split `Variables.f90` into focused domain/state modules.
  - [ ] Initialize scalar components deterministically.
  - [ ] Make modules `private` by default and expose explicit APIs.
  - [ ] Use `use ..., only:` imports.
- [ ] Reorganize `MOD_Recurrent_Functions.f90` into focused modules.
  - [ ] Normalization.
  - [ ] Geometry predicates/distances.
  - [ ] Results formatting.
  - [ ] Coarse mapping/setup.
  - [ ] Homogenization.
  - [ ] Modulation.
  - [ ] Fine-to-coarse current transfer.
  - [ ] Relaxation/SOR.
  - [ ] Frozen legacy coarse-input facade.
- [ ] Reformulate `Accel_Ini` behind a compatibility wrapper.
  - [ ] Separate frozen coarse-card parsing from validation and allocation.
  - [ ] Separate coarse topology/geometry construction.
  - [ ] Separate fine-to-coarse spatial and energy mapping.
  - [ ] Separate current-interface mapping and workspace allocation.
  - [ ] Redesign and verify `face_ID` and `face_ID_GL` assignment.
  - [ ] Define invariants for local/global face mapping and validate before indexing.
- [ ] Review and correct `dl_gl_Auto`, `dl_gl_Calc`, and related distance routines
  only after geometry tests freeze their intended behavior.
- [ ] Move SOR configuration/validation and the algorithm into a focused module;
  keep the arithmetic kernel simple/pure.

### P1.4 Consolidate iteration and strategy boundaries

- [ ] Introduce `linear_solver_options_t`, `linear_solver_result_t`, and reusable
  solver workspaces.
- [ ] Introduce concrete band and arrow operators with tested `apply` operations.
- [ ] Share one fine-mesh iteration controller between Quad and TriIso.
  - [ ] Keep topology-specific stencil assembly and current reconstruction concrete.
  - [ ] Share source, convergence, logging-status, and accelerator-dispatch logic.
- [ ] Share one coarse group/power iteration controller among acceleration variants.
- [ ] Replace parallel vector/matrix acceleration state with one authoritative state.
- [ ] Introduce an acceleration strategy selected once at orchestration:
  none, CMFD, pCMFD, HGCMFD, or pHGCMFD.
- [ ] Keep physics-specific coefficient/migration rules in concrete policies and
  verify them with conservation/parity tests.
- [ ] Reduce `Header.f90` to a thin composition root after ownership is established.

### P1.5 Complete output-I/O lifecycle refactoring

- [ ] Evolve `output_files_t` into an owned/injected output manager with output
  root, enable flags, open state, and error policy.
- [ ] Remove numerical-kernel dependence on the global `files` singleton.
- [ ] Make every `write(files%..., ...)` conditional on a valid enabled/open stream,
  or expose safe manager write methods.
- [ ] Prevent a failed input run from truncating a previous successful result.
- [ ] Make double-open, close, and disabled-output behavior deterministic.
- [ ] Do not create empty/unwritten output files unless they are intentional API.
- [ ] If errors must be persisted, add/manage an explicit error-file member;
  `error_unit` is process standard error, not a disk file by itself.
- [ ] Add result-output contract tests for dimensions, ordering, universes, flux,
  eigenfunction, power, and zero-power behavior.

### P1.6 Validate model invariants

- [ ] Require exactly one valid background figure and a valid material/universe for
  every physical cell.
- [ ] Require reduced IDs to be exactly `1..n_red` and outside IDs to be the defined
  sentinel.
- [ ] Validate positive volumes, face areas, distances, reciprocal neighbors, and
  opposite-face consistency.
- [ ] Validate fine-to-coarse spatial/energy mappings and HG background-node/arrow
  topology.
- [ ] Confirm scattering source/destination orientation before changing XS layout.
- [ ] Confirm whether figure-overlap first-match behavior is intentional.
- [ ] Correct the pi constant after a unit test; it is currently pi/16.
- [ ] Recheck and adjust Quad and TriIso convergence, boundary-condition,
  neighbor-face, and output handling with comparable sanity cases.

## P2 — Follow-up numerical and maintainability work

- [ ] Qualify `MOD_EigenvaluesSpectrumSolver.f90` separately.
  - [ ] Add 1x1, diagonal, symmetric/nonsymmetric 2x2, complex-pair, Arnoldi
    breakdown, QR factorization, and eigen-residual tests.
  - [ ] Fix undefined max-iteration status, NaN loop exit, work-vector shape,
    early-breakdown shapes, Ritz residual, dimension-one QR, and complex shifts only
    if the module will be retained.
  - [ ] Compare with a maintained LAPACK/Arnoldi implementation and decide whether
    to replace or remove the custom module from production.
- [ ] Move migration-matrix diagnostic formatting into a dedicated debug module.
  - [ ] Keep construction independent of formatting and file output.
  - [ ] Route output through the output manager's debug stream.
- [ ] Refactor `MOD_DebugChecks.f90`.
  - [ ] Separate reaction-rate calculation from diagnostic writing.
  - [ ] Remove full-size temporary arrays where streaming/small workspace suffices.
  - [ ] Make scattering headers group-count aware.
- [ ] Confirm pHGCMFD correction formulas, including the unused diffusion argument
  and hard-coded surface-flux term, with domain reference data.
- [ ] Decide and document whether `Crfc` `A_gl` means area or circumference.
- [ ] Document face-ID conventions for triangle/rhomboid/hexagon/octagon coarse
  shapes.
- [ ] Refactor Python utilities into import-safe libraries plus explicit CLI/config
  entry points after validating their indexing and output assumptions.
- [ ] Add professional project documentation: supported modes, build/test commands,
  numerical tolerances, failure statuses, input/output schemas, and architecture.

## P3 — Deferred input-reader modernization (not in first refactoring)

The following work is preserved from the previous checklist but is explicitly
deferred. Do not perform it while the input freeze is active.

### P3.1 Use managed units for input files

- [ ] In `MOD_Input_Reader.f90`, replace hard-coded unit 120 for
  `Options_Input.inp` with a local `options_unit` obtained through the chosen file
  manager.
  - [ ] Replace the open at the current line 30.
  - [ ] Replace every `read(120, ...)` and `rewind(120)`.
  - [ ] Replace early closes at the current lines 172 and 180 and final close at
    line 322 with one managed close path.
- [ ] Replace hard-coded unit 100 for `Material_Input.inp` with a local
  `material_unit`.
  - [ ] Replace the open at the current line 394.
  - [ ] Replace every `read(100, ...)` and `rewind(100)`.
  - [ ] Replace the close at the current line 526.
- [ ] Replace hard-coded unit 110 for `Geometry_Input.inp` with a local
  `geometry_unit`.
  - [ ] Replace the open at the current line 581.
  - [ ] Replace every `read(110, ...)` and `rewind(110)`.
  - [ ] Replace the close at the current line 1151.
- [ ] In `MOD_Recurrent_Functions.f90`, replace hard-coded unit 130 for both
  coarse-grid input paths with a local `grid_unit`.
  - [ ] CMFD: replace the open at the current line 286, every `read(130, ...)`, and
    close at line 470.
  - [ ] HGCMFD: replace the open at the current line 498, every `read(130, ...)`,
    rewind at line 530, and close at line 717.
- [ ] Define and document input-open/read failure propagation to callers.

### P3.2 Separate and harden parsing

- [ ] Separate mesh-generation routines from the input-card reader module.
- [ ] Check `iostat` and `iomsg` consistently for opens and reads.
- [ ] Handle malformed cards, premature EOF, missing sections, duplicate fields,
  invalid values, missing/duplicate background figures, and unmatched materials
  deterministically.
- [ ] Ensure every exit path closes its unit exactly once.
- [ ] Validate material-card group count against options.
- [ ] Replace whole-vector RM/TR/D absence heuristics with an explicit documented
  presence/default policy, with compatibility tests.
- [ ] Version any input syntax change and retain a documented compatibility path.

## P4 — Repository hygiene and retirement decisions

- [ ] Decide whether CMake is authoritative and whether the stale `.vfproj`, `.u2d`,
  and Intel wizard `ReadMe.txt` can be retired.
- [ ] Classify `HGCMFD/` as historical reference or migrate selected tested kernels.
  - [ ] Do not compile its byte-identical `MOD_HGCMFD_Package.f90` copy.
  - [ ] Do not adopt the decomposed header's inconsistent matrix/vector state flow.
- [ ] Decide whether empty `MOD_TriEq_Mesh.f90` and `MOD_lpCMFD_Package.f90` remain
  explicit unsupported placeholders or are removed until implementation.
- [ ] Define ignore/retention policy for `build/`, `x64/`, `.mod`, `.obj`, binaries,
  `__pycache__`, generated output, plots, and logs.
- [ ] Resolve the leading-space duplicate ` Coordinates_Vert.out` and other stale
  outputs only after provenance/retention decisions are made.
- [ ] Replace `builder.sh` with a fail-fast, documented build/test entry point or
  make it fail fast and call the canonical presets.

## Completed and verified in the current source

- [x] In `MOD_Quad_Mesh.f90`, replace the remaining `write(2, ...)` with
  `write(files%convergence, ...)`.
- [x] Correct `Results_Writing` to use its `n_tot` and `univ` dummy arguments rather
  than out-of-scope `Opt_lo`/`LO_Mesh` state.
- [x] Write the universe data once instead of writing the complete array `n_tot`
  times.
- [x] In `Header.f90`, pass `LO_Mesh(:)%univ` through an array constructor to
  preserve the universe output.

Completed implementation still requires the P1 output-contract tests before it is
treated as refactoring-safe behavior.
