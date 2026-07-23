# Implementation refinements

This document complements `todo.md`. It records distinct weaknesses in the
current solver implementation, their consequences, and possible remedies. It is
not an ordered work queue and does not replace the project roadmap.

Review snapshot:

- File: `src/solvers/MOD_BiCGSTAB_Solver.f90`
- SHA-256: `0652f416f5ca19871491da7d41a0604aa81546c2c8295970249a3333f67a096e`
- Scope: BiCGSTAB iteration, operator API, workspace, and preconditioners.
- Exclusion: the tri-, penta-, hepta-, EqTri-, and TriIso-style operator layouts
  are intentional designs and are not treated as problems here.

The fixed right-preconditioned BiCGSTAB recurrence has no identified sign or
operation-order error. The required refinements concern numerical scaling,
failure handling, stopping behavior, and API contracts.

Unless stated otherwise, verification cases should be added to
`tests/00_unit/TEST_BiCGSTAB_hybrid_suite.f90`. A successful solve must have a
finite solution and satisfy the stopping criterion using the explicitly
recomputed residual `b - A*x`.

## `MOD_BiCGSTAB_Solver.f90`

### R1 — Make the Euclidean norm overflow- and underflow-safe

- **Problem:** `euclidean_norm` first evaluates `dot_product(vector, vector)`.
  Squaring large components can overflow, while squaring small components can
  underflow to zero. The subsequent `max(0, ...)` does not repair either case
  and can hide a NaN on some compilers.
- **Impact:** The solver can report false convergence for a nonzero residual or
  convert an otherwise representable norm into a nonfinite value.
- **Possible solution:** Check that every component is finite, handle the zero
  vector separately, and use a scaled sum-of-squares calculation. One simple
  form is `scale*sqrt(sum((vector/scale)**2))`, where
  `scale = maxval(abs(vector))`. A LAPACK `xLASSQ`-style implementation avoids a
  second pass and provides the same protection. Do not assume that a compiler's
  `norm2` intrinsic uses scaled arithmetic without verifying it.
- **Verification:** Test zero, ordinary, NaN, infinity, values near `tiny()`, and
  values near `huge()`. Include vectors such as `[1e-200, 1e-200]` and
  `[1e200, 1e200]` when `prec` is binary64.

### R2 — Make breakdown detection scale-independent

- **Problem:** The rho, alpha-denominator, and `dot(t,t)` checks mix absolute
  comparisons with `tiny()` and one relative check. The product of two norms in
  `scaled_dot_is_zero` can itself overflow or underflow. The raw value of
  `omega` is compared with a dimensionless tolerance even though `omega` changes
  inversely when the matrix is rescaled. Only rho breakdown receives a restart.
- **Impact:** Multiplying an otherwise unchanged, well-conditioned system by a
  large or small constant can change convergence into a false breakdown.
  Certain nonsingular systems can also encounter a genuine BiCGSTAB shadow-vector
  breakdown with no recovery path.
- **Possible solution:** Express every dot-product breakdown test relative to
  the norms of its operands and evaluate the comparison with scaled arithmetic.
  Test the angle represented by `dot(t,s)` rather than the magnitude of `omega`,
  and obtain the `t` magnitude from the stable norm routine. On a recoverable
  alpha or omega breakdown, recompute the true residual and perform a bounded
  restart with a different deterministic shadow vector. If the restart repeats,
  return the specific breakdown status or use a documented fallback solver.
- **Verification:** Solve the same small system after scaling `A` and `b` over a
  wide exponent range; the status and final relative residual should remain
  consistent. Add explicit rho, alpha, and omega breakdown cases, including a
  nonsingular skew-symmetric case.

### R3 — Reject nonfinite component outputs immediately

- **Problem:** An operator or preconditioner can return success while placing
  NaN or infinity in its output. The solver trusts that output. Several derived
  scalars and updated vectors are also used before their finiteness is checked.
- **Impact:** A component defect can corrupt the recurrence or the solution and
  later be misreported as rho, alpha, or omega breakdown instead of a nonfinite
  value.
- **Possible solution:** After every successful operator or preconditioner call,
  verify that the entire output is finite before using it. Check `alpha`, `beta`,
  `omega`, and the solution update before continuing. Return
  `BICGSTAB_NONFINITE_VALUE`; reserve operator and preconditioner failure statuses
  for components that explicitly report failure. Native operator bind/setup
  routines should also reject nonfinite coefficients where practical.
- **Verification:** Use custom components that return NaN or infinity with a
  success status. Confirm immediate `BICGSTAB_NONFINITE_VALUE`, an unchanged
  solution when failure precedes an update, and no false breakdown status.

### R4 — Define and validate stopping-tolerance semantics

- **Problem:** The option checks reject negative tolerances but accept positive
  infinity. `absolute_tolerance + relative_tolerance*rhs_norm` is not checked for
  overflow. With the default zero absolute tolerance, a zero right-hand side
  requires an exactly zero residual, and `relative_residual` contains an absolute
  norm when `rhs_norm` is very small.
- **Impact:** An infinite tolerance can cause immediate false success. Zero- and
  small-right-hand-side behavior is surprising, and the reported relative
  residual changes meaning depending on the input scale.
- **Possible solution:** Require every real option and the computed residual limit
  to be finite and nonnegative. Document one reference scale for relative
  convergence. Either retain `||b||` and require callers to provide a useful
  absolute tolerance for zero `b`, or use a documented scale such as
  `max(||b||, ||r0||)`. Define `relative_residual` consistently; for a zero
  reference norm, report a documented sentinel or a separately named absolute
  measure.
- **Verification:** Reject NaN, infinity, and negative options. Test overflow in
  the residual-limit calculation and cover zero, subnormal, and ordinary
  right-hand sides with zero and nonzero initial guesses.

### R5 — Periodically replace the recursive residual

- **Problem:** The true residual is recomputed only after apparent convergence or
  at the iteration limit. During a long solve, the recursively updated residual
  can drift away from `b - A*x` because of roundoff.
- **Impact:** The solver can use inaccurate search directions, stagnate, waste
  iterations, or reach the iteration limit even though the true residual behaves
  differently.
- **Possible solution:** Add a configurable residual-replacement interval and an
  optional residual-gap threshold. Recompute `b - A*x` at that point; when the
  gap is significant, replace the recursive residual and restart the recurrence.
  Keep the existing true-residual verification before every success result.
- **Verification:** Use ill-conditioned and nonnormal systems that need many
  iterations. Record recursive and true residuals and confirm that replacement
  prevents an unbounded residual gap without changing easy-case convergence.

### R6 — Use local, finite-safe preconditioner pivot tests

- **Problem:** The default pivot threshold uses one global value proportional to
  `maxval(abs(diagonal))`. One very large diagonal entry can therefore reject a
  different, usable pivot. A zero override can permit reciprocal overflow, and
  ILU(0) has no shift or fallback when it encounters an unusable pivot.
- **Impact:** Jacobi, SGS, or ILU(0) setup can reject a safely solvable scaled
  system, or accept a pivot that generates infinity during setup or application.
- **Possible solution:** Validate that coefficients and the optional pivot
  tolerance are finite and that the tolerance is nonnegative. Compare each pivot
  with a row-local scale, then verify every reciprocal and computed factor is
  finite. For ILU(0), optionally support a documented diagonal shift, reordering,
  or fallback to a simpler preconditioner; otherwise return a specific pivot
  failure status.
- **Verification:** Cover diagonals with strongly different magnitudes, the
  smallest representable values, a nonsingular matrix with a zero leading ILU
  pivot, nonfinite coefficients, and every supported preconditioner.

### R7 — Make allocation failure recoverable and accurately reported

- **Problem:** Workspace allocations use `stat=`, but failure is reported as
  invalid input. Jacobi and ILU(0) setup allocations omit `stat=` and may
  terminate execution. An ILU(0) pivot failure can leave partial factor storage
  allocated.
- **Impact:** Resource exhaustion has inconsistent behavior and can terminate a
  library caller instead of returning a solver result.
- **Possible solution:** Check every allocation with `stat=` and optionally
  `errmsg=`, clean all partially initialized state on every failure path, and add
  a distinct allocation/resource status. Keep each preconditioner unready after
  any failed setup.
- **Verification:** Exercise failure paths with an injectable allocator or a test
  allocation limit. Confirm a deterministic status, no partially ready object,
  and successful setup when the same object is reused afterward.

### R8 — Make retained object lifetimes explicit and enforceable

- **Problem:** Banded operators retain pointers to coefficient dummy arguments,
  SGS retains an operator view, and user preconditioners retain procedure and
  context pointers. The current binder interfaces do not ensure that actual
  targets outlive those stored associations. Expressions, temporary contiguous
  copies, ordinary non-`TARGET` arrays, expired local contexts, and expired
  internal procedures can therefore leave dangling pointers.
- **Impact:** A valid-looking operator or preconditioner can later read undefined
  memory. The failure may depend on optimization level or call-site layout.
- **Possible solution:** Prefer owned allocatable coefficient storage. If zero-copy
  views are required, expose a clearly named view-binding API that accepts pointer
  actual arguments and documents the lifetime requirement; provide a separate
  copy/owning API for ordinary arrays. Apply an equivalent ownership rule to the
  callback context, or pass callback state directly to each solve instead of
  retaining it.
- **Verification:** Construct objects from local arrays, array sections,
  expressions, and long-lived pointer targets. After the construction scope ends,
  only owned objects or explicitly valid views should remain usable.

### R9 — Validate a custom operator before querying its order

- **Problem:** `bicgstab_core` calls `operator%order()` before
  `operator%is_valid()`. A custom invalid operator may not be able to compute its
  order safely.
- **Impact:** The solver can fault inside user code before returning
  `BICGSTAB_INVALID_INPUT`.
- **Possible solution:** Call `is_valid()` in a separate statement and return on
  failure before calling `order()`. Do not rely on logical-expression
  short-circuit evaluation, which Fortran does not guarantee. Also document that
  both queries must be side-effect-free for valid operators.
- **Verification:** Use an invalid custom operator with an order-call counter.
  Confirm that `order()` is not called and the solver returns invalid input.

### R10 — Enforce the fixed-preconditioner contract

- **Problem:** Standard right-preconditioned BiCGSTAB assumes one fixed linear
  preconditioning map. The public callback receives mutable context and can
  silently change its mapping between applications.
- **Impact:** A state-varying or nonlinear callback can invalidate the recurrence,
  producing unexplained stagnation or an incorrect convergence history even when
  every callback returns success.
- **Possible solution:** Document that mutation is allowed only for caches and
  diagnostics and that every application must represent the same linear map.
  Keep variable or nonlinear preconditioning outside this solver; if it is a
  requirement, provide a separate flexible Krylov method with an explicit API.
- **Verification:** Add a fixed user-preconditioner parity case. Add a deliberately
  changing callback as a contract test or documented unsupported example, rather
  than treating its behavior as valid BiCGSTAB output.
