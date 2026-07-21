# Refactoring preparation record

This document is the read-only analysis record for the planned first refactoring of
`Deterministic_FVM_Code`. It records the current architecture, numerical risks,
test strategy, and proposed modular/OOP boundaries so that the implementation can
be performed incrementally without repeatedly rediscovering the code.

## 1. Scope and constraints

### Current assignment

- Analyze the complete repository and the existing `todo.md`.
- Prepare, but do not perform, the first whole-code refactoring.
- Prefer OOP where it improves ownership, lifecycle management, or runtime strategy
  selection. Keep small stateless numerical kernels procedural.
- Improve readability, modularity, professional structure, and future extensibility.
- Give linear-solver correctness, tests, and measured optimization the highest
  priority.
- Plan a separate `test/` directory using the same TestDrive collector/driver idea
  as the HGCMFD reference testbed.
- Do not change the present input-reading behavior during the first refactoring.
- During this preparation pass, only `todo.md` may be changed among existing files;
  this new record is the only additional file.

### Interpretation of the input freeze

The frozen input subsystem is larger than `MOD_Input_Reader.f90`:

- `Options_Card`, `Material_Card`, and `Geometry_Card` are frozen.
- `Quad_Mesher`, `TriIso_Mesher`, the triangle construction helpers, and the
  disconnected `TriEq_Mesher` currently reside in the input module. They may be
  moved later without changing behavior, but their parsing-facing interfaces
  should initially remain compatible.
- The parsing section of `Accel_Ini` in `MOD_Recurrent_Functions.f90` is also input
  code. It reads `CMFD_Grid.inp` and `HGCMFD_Grid.inp` and must be behaviorally
  frozen in the first refactoring.
- A compatibility adapter may call the legacy readers and assemble their outputs
  into new aggregate objects. The card syntax, defaults, relative paths, and
  interpretation must not be redesigned in the first pass.

### Refactoring rule

Structural movement and numerical behavior changes must not be mixed in one step.
First characterize behavior; then fix one defect; then refactor behind passing
tests. Numerical equivalence, residuals, conservation identities, and explicit
failure statuses are the gates. Compilation alone is not a correctness gate.

## 2. Analysis status and evidence

- All active root Fortran source files were read and mapped.
- The dormant `HGCMFD/` decomposition was read and compared with the active root
  implementation.
- CMake, presets, shell build script, Intel Visual Studio metadata, Python tools,
  root input decks, all seven case directories, generated outputs, and PBR toy
  assets were inventoried and classified.
- The current `todo.md` was read in full and its completed and outstanding items
  were checked against the source.
- The reference testbed was inspected at
  `/home/guest/Davide/02_iMC/v_0.010/01_TESTBED_HGCMFD-implementation`.
  The actual directory uses lowercase `implementation`; the originally supplied
  uppercase spelling does not resolve on this case-sensitive filesystem.
- A compiler-only GNU Fortran audit was performed with module and temporary output
  directed to `/tmp`. It completed successfully. It reported many unused legacy
  declarations plus higher-value warnings about array temporaries, a guarded
  TriIso boundary index, real/complex conversion in the eigen code, and mixed-kind
  conversions. No repository output was produced by this audit.
- No solver or full application run was used to establish a baseline in this pass.
  Program startup opens output files with `status='replace'`, so an uncontrolled
  run could overwrite the existing output corpus.
- The local `.git` directory is empty and the directory is not a functional Git
  worktree. Existing binaries and output files therefore have no trustworthy
  source revision provenance.

## 3. Executive conclusions

1. Tests must precede the OOP reorganization. The current numerical code contains
   correctness defects that could otherwise be preserved or obscured by structural
   changes.
2. The highest-priority numerical defects are:
   - all three banded BiCGSTAB variants can turn an already converged solution into
     NaN;
   - Quad and TriIso inner iteration counters are never incremented, so their caps
     are ineffective;
   - dense Gaussian elimination has no singularity status and is repeatedly
     factorized for unchanged coarse matrices;
   - `Accel_Ini` declares `Opt_gl` as `intent(out)` and then reads its incoming
     values, which is undefined Fortran behavior;
   - the HGCMFD missing-grid fallback leaves identifiers undefined;
   - unsupported TriEq and lpCMFD choices are accepted as if operational.
3. OOP is valuable at ownership and strategy boundaries: simulation case, mesh,
   output lifecycle, linear solver, fine-mesh solver, and acceleration strategy.
   It is not valuable for every cell or every scalar formula.
4. The current per-cell allocatable layout already has poor locality. A hierarchy
   of polymorphic cell objects would make it worse. Aggregate objects should own
   contiguous arrays.
5. The active coarse solvers build both structured and dense representations, but
   use repeated dense Gaussian elimination. The HG arrowhead systems already have
   an inactive O(n) Schur solve. Correctness tests and benchmarks should determine
   the production backend.
6. `Service_Fcns` and `Accel_Ini` are the largest modularization seams, but their
   input-reading portions must initially remain behind a compatibility facade.
7. The split `HGCMFD/` tree is a useful responsibility map, not replacement code.
   It is unbuilt, type-incompatible with the root implementation, and contains a
   matrix/vector state mismatch.

## 4. Repository census

### 4.1 Active root Fortran

| File | Current responsibility | Refactoring status |
|---|---|---|
| `Header.f90` | Program composition, I/O startup, all input calls, mesh dispatch, result writing | Reduce to a thin composition root after aggregate ownership exists |
| `Variables.f90` | Kinds and all input, mesh, XS, coarse, and solver data types | Split into focused domain/state modules; default private |
| `MOD_IO_routines.f90` | Output file object plus global singleton | Evolve into an owned/injected output manager |
| `MOD_Input_Reader.f90` | Options/material/geometry parsing and fine mesh generation | Behavior frozen initially; later separate parsing, validation, and meshing |
| `MOD_Recurrent_Functions.f90` | Normalization, output, geometry, coarse input/setup, homogenization, modulation, currents, SOR | Split by responsibility; retain frozen reader facade |
| `MOD_BiCGSTAB_Solvers.f90` | Three duplicated structured BiCGSTAB solvers/matvecs | P0 characterization; consolidate after defects are fixed |
| `MOD_GaussElimination.f90` | Dense Gaussian elimination and pivoting | P0 tests/status; factor-once benchmark; likely reference backend |
| `MOD_EigenvaluesSpectrumSolver.f90` | Custom Arnoldi, matrix BiCGSTAB, complex QR | Compiled but unused; quarantine and test before retaining |
| `MOD_Quad_Mesh.f90` | Quad FVM assembly, source iteration, currents | Share iteration controller and pure kernels; retain topology-specific assembly |
| `MOD_TriIso_Mesh.f90` | TriIso FVM assembly, source iteration, currents | Same target as Quad; retain topology-specific assembly/current mapping |
| `MOD_CMFD_Package.f90` | CMFD coefficients, matrix assembly, power iteration | Acceleration strategy with shared controller |
| `MOD_pCMFD_Package.f90` | Partial-current CMFD variant | Acceleration strategy; physics-specific coefficient policy |
| `MOD_HGCMFD_Package.f90` | Heterogeneous global CMFD, dense and inactive arrow solvers | Unify state/operator; validate Schur path |
| `MOD_pHGCMFD_Package.f90` | Partial-current heterogeneous global CMFD | Same target as HGCMFD |
| `MOD_DebugChecks.f90` | Reaction-rate calculation and diagnostic writer | Compiled but uncalled; separate calculation from formatting/output |
| `MOD_TriEq_Mesh.f90` | Empty placeholder | Unsupported; reject explicitly until implemented |
| `MOD_lpCMFD_Package.f90` | Empty placeholder | Unsupported; reject explicitly until implemented |

The active/dormant Fortran and Python corpus is approximately 10,365 lines. The
largest active files are `MOD_Input_Reader.f90` (1,782 lines),
`MOD_Recurrent_Functions.f90` (1,344), and `MOD_TriIso_Mesh.f90` (886).

### 4.2 Dormant `HGCMFD/` tree

The directory has 21 Fortran files and approximately 1,692 lines. It is not added
to the root CMake build. All files use CRLF endings.

- `HGCMFD/MOD_HGCMFD_Package.f90` is byte-identical to the active root package.
  Recursive source discovery would create duplicate module/routine definitions.
- `MOD_HGCMFD_00-Header.f90` orchestrates the decomposed alternative.
- `MOD_HGCMFD_Variables.f90` defines separate `_2` data types that are not compatible
  with the root types.
- Coefficient/coupling modules:
  `MOD_HGCMFD_J-Homo-Cond.f90`, `MOD_HGCMFD_XS-Homo-Cond.f90`,
  `MOD_HGCMFD_J-to-El.f90`, `MOD_HGCMFD_D-tilde.f90`,
  `MOD_HGCMFD_D-hat.f90`, and `MOD_HGCMFD_MigrMatr.f90`.
- Iteration modules:
  `MOD_HGCMFD_S-fiss.f90`, `MOD_HGCMFD_S-tot.f90`,
  `MOD_HGCMFD_Update.f90`, `MOD_HGCMFD_Discrepancies.f90`,
  `MOD_HGCMFD_Inner-Iteration.f90`, and
  `MOD_HGCMFD_Outer-Iteration.f90`.
- Solver/support modules:
  `MOD_HGCMFD_Gauss-Elimination.f90`, `MOD_HGCMFD_BiCGSTAB.f90`,
  `MOD_HGCMFD_Shur-Complement.f90`, `MOD_HGCMFD_Normalization.f90`,
  `MOD_HGCMFD_Modulation.f90`, and `MOD_HGCMFD_SOR.f90`.

The decomposed header solves the matrix state `Serv_Matr` but then normalizes,
checks, and modulates `Serv_Vect%Phi`. This is a confirmed inconsistent-state
defect. The decomposition must not be copied unchanged.

### 4.3 Build and project metadata

- `CMakeLists.txt` builds one monolithic executable. There is no core library,
  CTest integration, warning profile, bounds checking, or floating-point trap
  profile. It creates `output/` in the source tree at configure time.
- `CMakePresets.json` provides Debug and Release binary directories/build types but
  does not pin a generator or compiler.
- `builder.sh` configures/builds Debug and Release but lacks fail-fast handling.
- `Deterministic_FVM_Code.vfproj` is stale Intel Visual Studio metadata. Debug x64
  selects ifx and Release x64 selects ifort. It omits `MOD_IO_routines.f90`, which
  the main program requires.
- `Deterministic_FVM_Code.u2d` is opaque Visual Studio conversion metadata, not
  source architecture.
- `ReadMe.txt` is Intel wizard boilerplate and names a nonexistent
  `Deterministic_FVM_Code.f90`.
- Existing `build/`, `x64/`, objects, modules, and executables are generated
  artifacts, not evidence of a reproducible clean build.

### 4.4 Python and auxiliary tools

- `Coordinates_Calculator.py` contains useful polygon/coordinate and plotting
  routines but executes a hard-coded case at import time, changes to a Windows
  path, and writes images. Later split reusable functions from a CLI/configuration
  entry point.
- `Trianglular_2D_plot.py` duplicates plotting helpers, uses a hard-coded Windows
  path, and executes file reads/plots at import time. Group slicing appears
  off-by-one, and one flux variable is overwritten with universe data. Treat this
  as low-priority tooling until validated.
- `20260211_PBR_Toy/Plotter.py` contains duplicated parsing for accelerated and
  non-accelerated logs, fixed-size assumptions, hard-coded Windows paths, and
  top-level plot generation. The `.png`, `.out`, and `Pebble_Inp.inp` files in that
  directory are experiment assets, not production FVM modules.

### 4.5 Input cases and generated data

The repository contains these case directories:

- `00_3x3_SMR_Reference`
- `00b_3x3_SMR_250x250`
- `01_3x3_Rhomboid_250x250`
- `02_3x3_Rhomboid_240x240`
- `02_SingleMaterial`
- `03_3x3_Octagon_250x250`
- `04_lpCMFD_SinglePin`

Together they cover Quad, TriIso/rhomboid, octagon, single-material, and
lpCMFD-labelled configurations. Several material and grid decks are byte-identical
between cases. These are useful fixtures, but not yet trusted numerical oracles.

The root `.out` files and `output/[Output]*.out` files contain conflicting result
sets. For example, root and `output/` eigenvalue histories differ materially and
have no recorded compiler/revision/deck provenance. They must not be declared
golden references without controlled reruns. Some output files are empty. There
are also two vertex files whose names differ only by a leading space:
`Coordinates_Vert.out` and ` Coordinates_Vert.out`.

Generated binaries, object/module files, plots, logs, and numerical outputs were
classified by origin and role. Binary/image content has no callable architecture;
it was not treated as source code.

## 5. Current execution and dependency flow

### 5.1 Program flow

`Header.f90` currently performs:

```text
open all output/debug files
→ Options_Card
→ Material_Card
→ Geometry_Card and fine-mesh construction
→ Accel_Ini when acceleration is enabled
→ FVM_Quad_Solver or FVM_TriIso_Solver
→ Results_Writing
→ close output/debug files
```

The TriEq case contains no solve. The program still proceeds to result writing,
where the flux is unallocated.

### 5.2 Dependency shape

```text
Variables
├── IO_module
├── Service_Fcns ───────────────┐
│   └── coarse input/setup      │
├── Input_Reader ───────────────┤
├── BiCGSTAB / Gauss / Eigen    │
├── CMFD-family accelerators ───┤
└── Quad/TriIso solvers ────────┤
                                ↓
                           Header program
```

Practical problems:

- Modules are public by default.
- Most `use` statements do not use `only` lists.
- Multiple accelerator modules export routines with identical names.
- `Input_Reader` imports all of `Service_Fcns` only to use a geometry predicate.
- Long-lived state is held as unrelated locals in `Header` and passed through large
  argument lists; there is no owner object.
- Quad/TriIso and the four accelerators duplicate iteration/source logic.

## 6. Current data model and contracts

### 6.1 Types

`Variables.f90` defines:

- `Int_input_check`, `Flo_input_check`, `Char_input_check`: value/presence wrappers.
- `Figure`: input geometry, vertices, material/universe, background marker.
- `Material`: multigroup cross sections.
- `XS_Data`: flattened cellwise cross-section fields.
- `Options_Data`: dimensions, flags, tolerances, iteration caps, and SOR settings.
- `LO_geom`, `LO_coeff`: fine geometry/topology/mapping and per-group coefficients.
- `GL_geom`, `GL_coeff`: coarse geometry/topology and coefficients.
- `Accel_Vars_Vect`, `Accel_Vars_Matr`: two representations of accelerator state.

Most scalar components have no default initialization. Missing or malformed input
can therefore propagate undefined identifiers.

The code mixes `dp`, `dp_var`, `doupr`, `real(8)`, and default-real literals. A
single `real64`-based kind module should eventually be authoritative.

`pi` is currently defined as `atan(1.d0)/4.d0`, which equals pi/16 rather than pi.
Its active effect is limited because its references are in dormant/incomplete
TriEq-related paths, but the constant is incorrect.

### 6.2 Dimensions and flattening

- Options input `nx` and `ny` are increased by two for a ghost layer.
- `Opt_lo%n_tot = n_x*n_y` includes ghost/outside cells.
- `Opt_lo%n_red` counts active cells inside the background figure.
- Fine arrays:
  - `size(LO_Mesh) = n_tot`
  - `size(LO_Coef) = n_tot*n_g`
  - `size(XS_lo%Tot) = n_tot*n_g`
- Intended group-major flattening:
  - full index: `(g-1)*n_tot + cell`
  - reduced index: `(g-1)*n_red + reduced_cell`
- `LO_Mesh%ID_Red = -1` identifies ghost/outside cells; positive IDs must be
  contiguous `1..n_red`.
- `LO_Mesh%lo_gl_Homo` maps fine cells to coarse regions; `-1` is the outside
  sentinel.
- `LO_Param((g-1)*n_tot+l)%En_ID` maps a fine group to a coarse group.
- `size(GL_Mesh) = Opt_gl%n_tot`; coarse coefficient/state arrays are group-major.
- `XS_Data%Scatt` is flattened on its first dimension and group-indexed on its
  second. Source/destination group orientation must be confirmed before a layout
  change.

These rules need named index helpers and invariant tests before reshaping arrays.

## 7. Input behavior to preserve initially

### 7.1 Options

`Options_Card` opens fixed relative `Options_Input.inp` on unit 120 and scans the
`% Local` and `% Acceleration` sections. It reads mesh type, dimensions, group
count, tolerances, iteration limits, SOR flag, and weight.

Current string-to-flag mapping:

| Input | Flag |
|---|---:|
| `Quad` | 0 |
| `Tri_Iso` | 1 |
| `Tri_Eq` | 2 |
| `None` | 0 |
| `CMFD` | 1 |
| `pCMFD` | 2 |
| `lpCMFD` | 3 |
| `HGCMFD` | 4 |
| `pHGCMFD` | 5 |

Parsing is case-sensitive, prefix/position dependent, and inconsistent in EOF and
error handling. These weaknesses are deferred, not to be silently corrected during
the first structural pass.

### 7.2 Materials

`Material_Card` opens fixed `Material_Input.inp` on unit 100. It allocates arrays
using the group count from options; the group count in the material deck is not
used as an independent authority.

Current whole-vector fallbacks:

- If the sum of removal XS is below `1e-12`, recompute all removal values.
- If the sum of transport XS is below `1e-12`, replace all values with total XS.
- If the sum of diffusion coefficients is below `1e-12`, compute
  `1/(3*transport)`.

Tests must characterize these exact heuristics before input modernization.

### 7.3 Geometry and fine mesh

`Geometry_Card` opens fixed `Geometry_Input.inp` on unit 110, parses figures and
boundary factors, allocates the fine arrays, invokes Quad or TriIso meshing, assigns
universe/material by center-point containment, identifies BC/current interfaces,
creates reduced IDs, expands material XS, and writes coordinates.

Implicit contracts include:

- one practical `!box` background;
- convex counter-clockwise polygons for the same-sign point test;
- first matching non-background figure wins when figures overlap;
- active neighbors are represented through generated ghost cells rather than
  arbitrary negative indices.

### 7.4 Coarse grids

`Accel_Ini` opens fixed unit 130:

- `CMFD_Grid.inp` for flags 1–3;
- `HGCMFD_Grid.inp` for flags 4–5.

It also allocates all coarse and accelerator storage, builds topology/geometry,
maps fine cells and energy groups, calculates faces, and initializes current maps.
This one routine therefore contains at least these future stages:

1. legacy coarse-card parse;
2. validation;
3. coarse allocation;
4. topology and geometry construction;
5. fine-to-coarse spatial mapping;
6. fine-to-coarse energy mapping;
7. current-interface mapping;
8. solver workspace creation.

Only stages 2–8 should be extracted in the initial pass if that can be done without
changing stage 1 behavior.

## 8. Output lifecycle

`output_files_t` owns unit numbers and type-bound open/close procedures, but not the
output path, enabled state, or error policy. A global singleton `files` is imported
throughout the code.

At startup, selected files are opened eagerly with `status='replace'`. This means a
later input failure can already have destroyed prior results. The object does not
remember which sets were enabled; close requires the caller to repeat booleans.
Calling open twice clears identifiers without first closing them. Disabling output
is unsafe because downstream code writes to units that remain `-1`.

The current output set includes log, convergence/eigenvalue history, flux,
eigenfunctions/eigenvalues, performance, coordinates, power, and universes. The
debug set includes XS, flux/current/coefficient histories, migration matrix, and
iteration diagnostics. Several streams are opened but never written.

Future `output_manager_t` responsibilities:

- own output root and enable flags;
- be safely openable/closable once;
- expose query/write methods or valid optional sinks;
- avoid truncating a successful result until input validation succeeds;
- separate calculation from formatting;
- propagate a structured status instead of mixing print-and-continue and
  `error stop`.

Input-unit migration belongs to the deferred input phase even if the same low-level
file helper is reused.

## 9. Numerical solver audit

### 9.1 Structured BiCGSTAB

`MOD_BiCGSTAB_Solvers.f90` has three near-identical implementations:

- `BiCGSTAB_Classic` / `matvect_penta`
- `BiCGSTAB_EqTri` / `matvect_penta_EqTri`
- `BiCGSTAB_TriIso` / `matvect_penta_TriIso`

Confirmed issues:

- no initial residual check;
- the `s`-residual early-success block is commented out;
- an exact initial guess or one-step solution can proceed to `0/0` in omega and
  corrupt a valid solution with NaN;
- no breakdown guard for rho, omega, `dot(r_hat,v)`, or `dot(t,t)`;
- convergence uses relative change in `x`, not the normalized residual;
- zero RHS and zero solution make relative solution denominators unsafe;
- status is only success/failure, with no iterations, residual, or reason;
- seven work arrays are allocated on every call;
- the three algorithm bodies are duplicated;
- explicit-shape matvec arguments trigger array temporaries in some calls;
- `matvect_penta` unconditionally addresses neighbors of the first and last cell,
  making very small grids unsafe even when a coefficient is zero;
- the guarded TriIso `x(i-1)` path receives a compiler bounds warning and needs a
  boundary fixture under runtime checking.

Target: one tested BiCGSTAB iteration core receiving an operator application
procedure or minimal operator object, reusable workspace, and a structured result.

### 9.2 Dense Gaussian elimination

`GE_Expl_Main` copies the full matrix/RHS, performs partial pivoting, and back
substitution.

Confirmed issues:

- no singular or near-singular pivot detection;
- no returned status;
- arithmetic row swaps can overflow/cancel;
- precision declarations are inconsistent;
- a full O(n^3) factorization is repeated inside each unchanged CMFD-family group
  sweep;
- input immutability and n=1 behavior are undocumented.

Tests should establish a dense reference solver first. Then benchmark:

- factor once / solve many using owned LU factors;
- LAPACK `GETRF/GETRS` or `GESV`, if an external dependency is accepted;
- structured band or arrow alternatives.

Dense solve should likely remain a reference/debug backend for arrow systems.

### 9.3 HG arrowhead solvers

HGCMFD/pHGCMFD build an arrowhead matrix per group: one background row, one
background column, and diagonal remainder. The source already contains inactive
Schur and arrow-BiCGSTAB routines.

- The Schur solve is O(n) and is the first optimized path to validate.
- It requires explicit checks for zero diagonal entries and a zero/near-zero Schur
  complement.
- Arrow BiCGSTAB has the same missing initial convergence/breakdown/status problems
  as the banded implementations.
- Dense and structured migration assembly must be checked for exact algebraic
  parity before changing the backend.

### 9.4 Eigenvalue-spectrum module

`MOD_EigenvaluesSpectrumSolver.f90` is compiled but has no active caller. It should
not be integrated into the refactored runtime until a separate test decision is
made.

Confirmed issues include:

- matrix BiCGSTAB leaves `info` undefined on maximum iterations and does not exit
  its NaN branch;
- Arnoldi allocates a work vector with reduced dimension instead of operator
  dimension, hidden because the current caller uses equal dimensions;
- early Arnoldi breakdown uses incompatible shapes;
- the Ritz residual uses the wrong attained dimension/component;
- a zero starting vector is normalized without a guard;
- QR for dimension one indexes element zero;
- the two-by-two shift omits the second off-diagonal of a nonsymmetric matrix;
- complex roots are assigned through real variables and imaginary parts are lost;
- accumulated QR vectors are Schur vectors, not general nonsymmetric eigenvectors;
- full dense arrays and potentially one million dense QR steps make the approach
  expensive.

For small dense matrices, benchmark a maintained LAPACK eigensolver. For large
problems, prefer a maintained Arnoldi implementation rather than repairing this as
part of the general refactor.

## 10. Fine-mesh and acceleration audit

### 10.1 Quad and TriIso

Both modules contain nearly the same outer/inner source-iteration controller,
source functions, flux extension, accelerator selection, logging, and failure
handling. Topology-specific matrix/current code differs and should remain concrete.

Confirmed issues:

- `iter_inn_lo` is initialized but never incremented in both inner loops;
- solver failure is not represented by a structured result;
- maximum outer iteration can still be reported as convergence;
- final extended flux can be stale if acceleration changes reduced flux on the last
  permitted outer iteration;
- matrix stencils rely on undocumented reduced-ID ordering;
- neighbor IDs are indexed before explicit validation;
- matrix assembly and current reconstruction apply boundary predicates
  inconsistently;
- same-named helper routines are public from both modules.

Target: one shared non-polymorphic iteration engine with topology-specific operator
assembly/current reconstruction strategies. Pure source and coefficient formulas
remain ordinary module procedures.

### 10.2 CMFD family

The four active acceleration packages share this broad lifecycle:

```text
fine-to-coarse currents
→ homogenize XS/flux
→ build D-tilde and correction coefficients
→ assemble structured and dense migration operators
→ coarse inner/outer iteration using repeated dense GE
→ normalize
→ positivity gate
→ modulate fine flux
→ optional SOR
```

Confirmed issues:

- all use `product(Phi)>0` as a positivity gate; this can underflow, overflow, and
  accepts an even number of negatives;
- correction formulas divide by flux, sums, areas, distances, or homogenized
  integrals without consistent guards;
- inner/outer failure status is not propagated;
- both dense and structured storage are built even though only dense GE is active;
- dense storage is allocated as `(n_coarse*n_group)^2` despite group-block solves;
- matrices are unchanged during a solve but repeatedly factorized;
- HG and pHG call D-hat construction twice with initially identical copied flux;
- background node 1 and arrow/star connectivity are implicit HG invariants;
- `pHGCMFD_D_Hat_Build` has an unused diffusion argument and a hard-coded zero
  surface flux; this requires physics confirmation, not an assumed cleanup.

Target: a common acceleration interface and shared power/group iteration controller,
with concrete coefficient/migration policies and an injected tested solver.

### 10.3 Shared services

`Service_Fcns` mixes unrelated concerns. Proposed procedural modules:

- `normalization_m`
- `geometry_predicates_m`
- `results_format_m`
- `coarse_mapping_m`
- `homogenization_m`
- `modulation_m`
- `current_transfer_m`
- `relaxation_m`
- a frozen `legacy_coarse_input_m` facade

Known guards/invariants to add after characterization:

- nonzero volume and normalization factor;
- nonzero homogenized flux integral;
- nonzero modulation denominator;
- successful `findloc` before indexing;
- positive face area/distance;
- valid fine/coarse map ranges;
- finite elementwise-positive flux where required;
- zero-power behavior in `Results_Writing`.

## 11. Confirmed defect register

| ID | Priority | Location | Defect/risk | Required gate before change |
|---|---|---|---|---|
| D-01 | P0 | `MOD_BiCGSTAB_Solvers.f90` | Exact/one-step convergence can become NaN; no breakdown guards | BiCG residual and edge-case suite |
| D-02 | P0 | Quad/TriIso inner loops | Counter never increments; cap ineffective | Iteration-limit test with timeout |
| D-03 | P0 | `MOD_GaussElimination.f90` | Singular pivots divide silently; unsafe row swap | Dense solver success/failure suite |
| D-04 | P0 | all active accelerators | Repeated O(n^3) factorization of unchanged matrices | Factor-count/time/residual benchmark |
| D-05 | P0 | all active accelerators | `product(Phi)>0` is not an elementwise finite positivity test | Flux validation tests |
| D-06 | P0 | `Accel_Ini` | `Opt_gl intent(out)` is read on entry | Coarse-initialization characterization |
| D-07 | P0 | HG fallback | Missing grid leaves `box_mat_ID`/related mapping undefined | Missing-grid integration test |
| D-08 | P0 | main/options | TriEq and lpCMFD accepted but unimplemented | Explicit rejection tests |
| D-09 | P0 | no-acceleration call | Unallocated `Elem_CMFD` passed to assumed-shape dummy | No-acceleration checked-build case |
| D-10 | P1 | iteration controllers | Failure/cap can be logged as convergence; stale final extension | Status/log/final-flux tests |
| D-11 | P1 | normalization/coupling | Multiple unchecked zero denominators and `findloc(…)=0` indices | Kernel failure/invariant suites |
| D-12 | P1 | `Variables.f90` | pi is pi/16; mutable, mixed precision | Constant/unit test; TriEq decision |
| D-13 | P1 | coarse geometry | HG reorder and face indexing assume background is first | Permuted-input mapping test |
| D-14 | P1 | `Results_Writing` | Relative power has zero-count/sum divisions | Zero-fission/output contract test |
| D-15 | P2 | eigen module | Multiple correctness, shape, complex, status, performance defects | Separate eigen qualification suite |
| D-16 | P2 | dormant HG split | Matrix state solved, vector state normalized/modulated | Do not adopt; parity tests before reuse |
| D-17 | P2 | debug checks | Excess full-size temporary storage and fixed two-scatter header | Debug output characterization |

Line numbers in source should be treated as snapshot references, not stable task
identifiers. Use routine names and defect IDs in future work.

## 12. TestDrive reference and new test architecture

### 12.1 Reference pattern inspected

The usable reference is the per-unit test configuration below
`01_TESTBED_HGCMFD-implementation/tests/00_unit/HGCMFD`, not the combined
`tests/CMakeLists.txt`.

The reference pattern is:

- a vendored `testdrive.F90` (reported TestDrive 0.5.0);
- one suite module per concern;
- each suite exposes a collector returning `unittest_type(:)` via
  `new_unittest`;
- one wrapper registers ordered `new_testsuite` values, calls `run_testsuite`,
  accumulates failures, and terminates nonzero;
- CMake compiles TestDrive with `WITH_QP=0` and `WITH_XDP=0`;
- a separate Fortran module output directory is used;
- the driver is registered with CTest and runs with the fixture directory as its
  working directory.

Reference coverage includes input conversion, D-tilde/D-hat, dense/vector migration
parity, one nominal arrow BiCG solve, one pivoted dense solve, inner iterations,
normalization, and one one-group reconstruction. It is useful precedent, but it is
not a complete solver qualification suite.

Do not copy these reference defects:

- stale incomplete combined `tests/CMakeLists.txt`;
- direct source globbing instead of linking a production library;
- hard-coded absolute/source-tree fixture paths;
- fixed I/O units that could race under parallel execution;
- one CTest entry with no labels/timeouts/granularity;
- API assumptions from a different experimental implementation;
- nominal-only solver tests;
- untracked/dirty build artifacts as proof of success.

### 12.2 Proposed `test/` layout

The requested directory name is singular:

```text
test/
├── CMakeLists.txt
├── testdrive/
│   └── testdrive.F90
├── driver/
│   └── main_test_driver.f90
├── unit/
│   ├── test_stencil_matvec_suite.f90
│   ├── test_bicgstab_suite.f90
│   ├── test_dense_solver_suite.f90
│   ├── test_arrow_solver_suite.f90
│   ├── test_normalization_suite.f90
│   ├── test_geometry_mapping_suite.f90
│   └── test_io_suite.f90
├── integration/
│   ├── test_iteration_suite.f90
│   ├── test_acceleration_suite.f90
│   └── test_sample_cases_suite.f90
├── fixtures/
│   ├── solver/
│   ├── mesh/
│   └── cases/
└── benchmark/
    ├── benchmark_linear_solvers.f90
    └── README.md
```

Implementation rules:

- Build a production `fvm_core` library and link both the executable and tests to
  it; do not compile a divergent copy of production sources in each test target.
- Use explicit source lists.
- Copy/configure fixtures into the test build tree.
- Isolate `.mod` files per target/configuration.
- Run I/O or end-to-end tests in isolated temporary working directories because
  current output startup replaces files.
- Disable TestDrive parallel execution for suites that touch legacy fixed units,
  until those units are removed.
- Add a strict Debug target with bounds checks, backtraces, warnings, floating-point
  traps, and compiler-specific equivalents.
- Keep performance measurements out of hard pass/fail CTest thresholds at first.

### 12.3 P0 solver test matrix

#### Operator application

- Compare Classic, EqTri, TriIso, and arrow matvec results against independently
  assembled dense matrices.
- Test arbitrary vectors and every basis vector.
- Cover first/last row, corners, edges, rectangular grids, every TriIso orientation,
  n=1/n=2 where meaningful, and fixed-seed valid random coefficients.
- Run with bounds checking.

#### BiCGSTAB

- exact initial solution;
- zero RHS/zero solution;
- one-step diagonal solution;
- known nonsymmetric diagonally dominant systems;
- nonzero initial guess;
- scaled matrices/RHS;
- max iterations zero/reached;
- singular/zero operator;
- explicit rho/alpha/omega/t-norm breakdown cases;
- nonfinite input.

Acceptance on success:

```text
all outputs finite
and
||b - A*x||_2 <= atol + rtol*||b||_2
```

Failure cases must return a deterministic status and reason without NaN/Inf or an
infinite loop.

#### Dense solver

- 1x1, identity, pivot-required, nonsymmetric well-conditioned, and fixed-seed
  diagonally dominant systems;
- exact singular and near-singular systems;
- large/small scale ranges;
- input `A`/`b` remain unchanged when promised;
- residual and known-solution checks;
- deterministic singular status.

#### Arrow/solver parity

- validate structured matvec against dense conversion;
- validate dense and vector migration assembly parity;
- compare Schur, dense reference, and arrow BiCGSTAB solutions/residuals;
- cover n=1, n=2, nominal medium, and production-like sizes;
- test zero arrow diagonals and near-zero Schur complement.

#### Eigen qualification, if retained

- 1x1, diagonal, symmetric 2x2, nonsymmetric real pair, and complex conjugate pair;
- early Arnoldi breakdown and zero start vector;
- unordered eigenvalue-set checks;
- `||A*v-lambda*v||` checks;
- QR factorization and unitary/orthogonal identity checks.

### 12.4 Integration and invariant tests

#### Iteration and error propagation

- inner limit reached exactly, without a hang;
- outer limit status;
- linear-solver failure propagation;
- success is not logged after failure/cap;
- returned extended flux matches final reduced flux;
- finite eigenvalue and convergence metrics.

#### Geometry/topology

- one background figure, unique materials/IDs;
- positive cell volumes, face areas, and distances;
- valid neighbor IDs and reciprocal opposite faces;
- reduced IDs exactly `1..n_red`;
- Quad physical volume sum equals box area;
- all four TriIso orientations and reciprocal adjacency;
- point-in-polygon inside/outside/edge behavior;
- valid fine-to-coarse spatial/energy mappings;
- every required `findloc` face match is positive;
- HG background-node/arrow topology invariant;
- current cancellation on paired interior faces.

#### XS and coupling

- expected shapes and finite values;
- frozen RM/TR/D fallback behavior;
- group-major material-to-cell expansion;
- homogenization preserves volume-integrated reaction rates;
- normalization postcondition and idempotence;
- zero denominator returns defined failure;
- modulation identity when coarse reference/current averages match;
- integrated coarse current equals contributing fine current integral;
- scattering source/destination convention, after domain confirmation.

#### I/O and output

- enabled units are valid after open and reset after close;
- double-open and mismatched close flags do not leak units;
- output-disabled mode never writes an invalid unit;
- defined behavior when output directory is absent;
- failed input does not destroy a previous successful result;
- result dimensions, universe count, flux ordering, and zero-power behavior;
- no empty files unless part of the intended interface.

#### Controlled sample baselines

Run every supplied case in a copied temporary case directory. Record:

- source snapshot/commit identifier;
- compiler, version, flags, and build type;
- SHA-256 of input decks;
- exit/status;
- final eigenvalue, iteration counts, residuals;
- flux/power norms and conservation metrics;
- selected numerical outputs with explicit tolerances.

Do not bless the present root/output `.out` corpus without this provenance.

## 13. Performance investigation plan

Correctness and performance results must be reported separately. Initial
benchmarks should not fail CI on wall-clock thresholds.

Measure representative small, medium, and production-like coarse problems:

- current repeated dense GE;
- factor-once/solve-many dense LU;
- LAPACK dense/banded LU if dependency policy allows;
- current structured BiCGSTAB;
- BiCGSTAB with reusable workspace and candidate preconditioner;
- HG arrow Schur solve;
- HG arrow BiCGSTAB.

Record:

- normalized residual/backward error;
- convergence/failure status and iterations;
- factorization count;
- median/min timings after warm-up;
- allocation count where measurable;
- peak memory/bytes by representation;
- scaling with cells/groups;
- compiler and flags.

Specific candidates:

- factor unchanged coarse matrices once;
- reuse the seven-vector BiCG workspace;
- stop constructing dense and structured matrices simultaneously unless requested;
- store one group block rather than a full multigroup dense matrix;
- exploit arrow O(n) storage/solve;
- replace per-cell allocatables with contiguous aggregate arrays in a later,
  separately benchmarked data-layout step;
- benchmark type-bound/abstract matvec dispatch against direct procedure calls
  before using dynamic dispatch inside hot loops.

## 14. Proposed target architecture

### 14.1 Design principles

- OOP for ownership, lifecycle, validation, and interchangeable strategies.
- Procedural/pure code for stateless formulas and tight kernels.
- One authoritative state representation per algorithm.
- Contiguous arrays inside aggregate types; no polymorphic object per cell.
- `private` by default and explicit small `public` APIs.
- `use ..., only:` imports.
- Dependency direction toward domain/numerics, never from low-level kernels to the
  program/global output singleton.
- Preserve existing array layouts during the first mechanical extraction unless a
  dedicated test/benchmark justifies changing them.

### 14.2 Proposed components

#### Foundation

- `kinds_m`: `real64` precision and immutable constants.
- `status_m`: status/error code and message conventions.
- `indexing_m`: full/reduced/group index helpers and checked conversions.

#### Domain ownership

- `simulation_case_t`
  - local/global options;
  - material library and XS fields;
  - fine mesh/model;
  - optional coarse model;
  - acceleration and linear-solver workspaces;
  - result/history object.
- `simulation_t`
  - `initialize_from_legacy_input`;
  - `validate`;
  - `run`;
  - `write_results`/`finalize`.
- `simulation_result_t`
  - flux, eigenvalue, power, histories, timings, final status.

#### Mesh/model

- `fine_mesh_t`/`coarse_mesh_t` own contiguous topology and geometry arrays and
  expose validation.
- Quad and TriIso concrete topology builders/assemblers.
- `cross_section_field_t` owns shaped data and supplies temporary legacy flattened
  views/adapters while interfaces migrate.

#### Linear algebra

- `linear_solver_options_t`: tolerances and limits.
- `linear_solver_result_t`: code, reason, iterations, residual.
- `bicgstab_workspace_t`: reusable vectors and `ensure_size`.
- `band_operator_t`: coefficients/topology and `apply`.
- `arrow_operator_t`: coefficients, validation, `apply`, optional `to_dense` for
  tests/debug.
- `dense_lu_t`: owned factors/pivots with `factorize` and `solve`.
- A minimal abstract `linear_operator_t`/`linear_solver_t` is acceptable only if
  multiple runtime strategies need one API and dispatch overhead is measured.

#### Physics iteration

- shared fine-mesh iteration controller;
- Quad/TriIso operator/current policies;
- shared coarse group/power iteration controller;
- one state object rather than parallel vector/matrix states.

#### Acceleration

- strongest justified polymorphic seam:
  - no acceleration;
  - CMFD;
  - pCMFD;
  - HGCMFD;
  - pHGCMFD.
- Select once through a factory at the composition boundary.
- Reject lpCMFD until implemented.
- Concrete variants retain physics-specific coefficient/migration rules while
  sharing controller, validation, solver, normalization, modulation, and status.

#### I/O

- `output_manager_t` owns path, flags, units, and lifecycle.
- Inject it into orchestration; numerical kernels return data/status and do not
  write files.
- Keep the legacy input reader as an adapter in the first pass.

### 14.3 Target dependency direction

```text
kinds/status/indexing
        ↓
domain arrays and mesh types
        ↓
pure geometry/source/coupling kernels
        ↓
operators and linear solvers
        ↓
fine iteration + acceleration strategies
        ↓
simulation/application service
        ↓
legacy input adapter + output manager
        ↓
thin main program
```

## 15. Recommended staged implementation

### Stage 0: safety net

1. Establish a usable source-control baseline or equivalent immutable snapshot.
2. Split CMake into `fvm_core` plus the existing executable without changing
   runtime behavior.
3. Add `test/`, vendored/pinned TestDrive, CTest, isolated module/output dirs, and
   a strict diagnostic build.
4. Add matvec, linear-solver, iteration-limit, and coarse-initialization tests.
5. Create controlled integration baselines from copied case directories.

### Stage 1: isolated correctness repairs

1. Fix BiCG initial convergence, `s` convergence, breakdown handling, residual
   criterion, and structured result.
2. Increment and correctly report fine inner/outer iteration status.
3. Add dense pivot checks/status and safe row exchange.
4. Correct `Accel_Ini` intent/undefined-state use and HG fallback behind tests.
5. Explicitly reject unsupported TriEq/lpCMFD.
6. Guard no-acceleration optional coarse arguments.
7. Replace product-based positivity with finite elementwise validation.
8. Add denominator/map guards and defined failure propagation.

Each repair should be isolated from structural moves and update baselines only when
the behavior change is intentional and documented.

### Stage 2: solver optimization based on measurements

1. Add factor-once/solve-many state.
2. Validate/enable arrow Schur for HG variants if residual and breakdown tests pass.
3. Decide dense/banded/iterative backends by problem shape and benchmark results.
4. Reuse workspaces and avoid constructing unused representations.

### Stage 3: aggregate ownership and modular extraction

1. Introduce `simulation_case_t` around current arrays.
2. Split foundational/domain types while preserving layout.
3. Extract pure kernels from `Service_Fcns`.
4. Split `Accel_Ini` stages behind the frozen parse wrapper.
5. Introduce solver result/options/workspace objects.
6. Consolidate Quad/TriIso controller logic.
7. Consolidate coarse controller logic and introduce accelerator strategies.
8. Make modules private by default and narrow imports.
9. Reduce `Header.f90` to orchestration.

### Stage 4: I/O and tooling

1. Replace global output singleton with an injected manager.
2. Separate diagnostics/calculation/formatting.
3. Clean build metadata, generated-artifact policy, and documentation.
4. Refactor Python utilities into import-safe libraries and CLI entry points.

### Deferred stage: input redesign

Only after the first refactor is stable:

- migrate hard-coded units 120, 100, 110, and 130;
- standardize `newunit`, `iostat`, `iomsg`, and close behavior;
- validate missing/duplicate/invalid fields deterministically;
- separate card parsing from mesh generation;
- replace integer flags/strings with validated enums/factories while maintaining a
  compatible input syntax or versioned migration path.

## 16. Refactoring acceptance gates

A structural step is acceptable only when:

- strict Debug compilation succeeds;
- unit/integration tests pass;
- no new warnings in changed code;
- solver successes meet residual tolerances;
- failures return defined statuses without NaN/Inf/hangs;
- controlled sample-case eigenvalue, flux, power, and conservation metrics remain
  within documented tolerances unless an intentional defect fix explains a change;
- performance-sensitive changes include benchmark evidence;
- input files and parsing behavior remain unchanged in the first refactor;
- source responsibility and public API are smaller and clearer than before.

## 17. Open decisions requiring confirmation

Do not guess these during implementation:

1. Is CMake the sole authoritative build, allowing the stale Intel `.vfproj` to be
   retired, or must both remain supported?
2. Is the split `HGCMFD/` tree historical/reference material, or must it be migrated
   into the production design?
3. What is the authoritative source-to-destination convention for `XS_SC`?
4. For the `Crfc` coarse figure, is `A_gl` intended to be area or circumference?
   The comment and formula disagree.
5. Are hard-coded face-ID conventions for triangle/rhomboid/hexagon/octagon inputs
   part of the external contract?
6. Is dense GE meant as the production coarse solver, or only as a reference while
   structured solvers become primary?
7. For overlapping non-background figures, is first-match input order intentional?
8. Must output formatting be byte-identical, or is numerical equivalence with a
   documented schema sufficient?
9. Should the unused custom eigen-spectrum module be retained, replaced with a
   maintained library, or removed from the production build?
10. What failure policy is required: returned status, top-level `error stop`, or
    recoverable per-case errors?

## 18. Routine reference map

This is the current routine-level navigation map. It is a snapshot; names are more
stable references than line numbers during the refactoring.

### `MOD_IO_routines.f90` / `IO_module`

- `open_file_main`: clears IDs, opens enabled output/debug sets, writes the banner.
- `file_ID_clear`: resets every stored unit to `-1` without closing an existing
  unit.
- `open_output_files_t`: opens the ten standard output files.
- `open_debug_files`: opens the thirteen debug files.
- `open_single_file_proc`: `newunit` wrapper that currently terminates on failure.
- `close_file_proc`: closes one unit after an `inquire` check.
- `write_log_header`: terminal/file banner.
- `close_file_main`: closes sets according to caller-supplied booleans.

### `MOD_Recurrent_Functions.f90` / `Service_Fcns`

- `Flux_Vol_Normalization`: global volume normalization.
- `Flux_lo_Vol_Normalization`: reduced fine-flux normalization; no active caller
  was found.
- `Output_Converter`: reduced-to-full flux expansion; active use is commented out.
- `Results_Writing`: universe, flux, dominant eigenfunction, raw power, and relative
  power output.
- `dl_gl_Auto`, `dl_gl_Calc`: point-to-edge distance helpers.
- `Point_in_Polygon_Test_2D`: same-sign edge test for convex CCW polygons.
- `Accel_Ini`: coarse parse/allocation/topology/mapping/workspace routine.
- `Distance`: Euclidean 2-D distance.
- `Homogenization_DerType`: fine-to-coarse XS and flux homogenization.
- `Modulation`: applies a coarse correction ratio to fine flux.
- `Currents_lo_to_gl_DerType`: integrates fine currents on coarse faces.
- `SOR_Accel`: linear relaxation update.

### `MOD_Input_Reader.f90` / `Input_Reader`

- `Options_Card`: option parsing, ghost-size adjustment, string-to-flag mapping.
- `Material_Card`: material/XS parsing and RM/TR/D fallbacks.
- `Geometry_Card`: geometry/BC parsing, fine allocation, mesher dispatch,
  material/BC/reduced mapping, XS expansion, coordinate output.
- `Quad_Mesher`: Cartesian mesh including a ghost layer.
- `TriIso_Mesher`: isosceles-triangle mesh.
- four triangle construction helpers used by `TriIso_Mesher`.
- `TriEq_Mesher`: disconnected equilateral-triangle prototype.

### `MOD_BiCGSTAB_Solvers.f90` / `BiCGSTAB_Solvers`

- `BiCGSTAB_Classic` and `matvect_penta`.
- `BiCGSTAB_EqTri` and `matvect_penta_EqTri`.
- `BiCGSTAB_TriIso` and `matvect_penta_TriIso`.

The three solver bodies differ primarily in the operator application. This is the
main duplication to remove after the common algorithm is qualified.

### `MOD_GaussElimination.f90` / `GaussElimination`

- `GE_Expl_Main`: matrix/RHS copy, elimination, back substitution.
- `partial_pivot`: pivot search and current arithmetic row swap.

### `MOD_EigenvaluesSpectrumSolver.f90`

- `Eigen_Routine`: public orchestration.
- `Arnoldi_alg`: basis/Hessenberg construction and reduced eigenproblem.
- `norm_l2`: vector norm helper.
- `BiCGSTAB_Mat`: dense-matrix BiCGSTAB used inside Arnoldi.
- `QR_main_Complex`, `QR_shift_Complex`, `QR_decomp_Complex`: custom complex QR
  iteration and factorization helpers.

### `MOD_Quad_Mesh.f90` / `FVM_Quad`

- `FVM_Quad_Solver`: outer/inner source iteration and accelerator dispatch.
- `FVM_Quad_MatrixBuilder`: structured operator coefficients.
- `Int_Coeff_D`, `Phi_Sup`: interface/boundary helper formulas.
- `FissionSourceDistribution`, `Total_Source_SingleGroup`, `Total_Source`: source
  kernels.
- `Local_Currents_Quad`: face-current reconstruction.
- `Phi_Extender`: reduced-to-full flux mapping.

### `MOD_TriIso_Mesh.f90` / `FVM_TriIso`

- `FVM_TriIso_Solver`: outer/inner source iteration and accelerator dispatch.
- `FVM_TriIso_MatrixBuilder`: orientation-dependent structured operator.
- duplicated coefficient/source helpers corresponding to the Quad routines.
- `Phi_Extender`: reduced-to-full flux mapping.
- `Local_Currents_TriIso`: triangle face-current reconstruction.

### Active acceleration packages

- `MOD_CMFD_Package.f90`: `CMFD_Header`, D-tilde/D-hat construction, implicit and
  explicit migration assembly.
- `MOD_pCMFD_Package.f90`: `pCMFD_Header`, D-tilde, positive/negative
  partial-current corrections, implicit and explicit migration assembly.
- `MOD_HGCMFD_Package.f90`: `HGCMFD_Header`, element current/coefficient
  condensation, star/arrow and dense migration assembly, inactive optimized-D
  experiments, Schur solve, arrow BiCGSTAB, and arrow matvec.
- `MOD_pHGCMFD_Package.f90`: `pHGCMFD_Header`, corresponding partial-current
  element/migration routines, Schur solve, arrow BiCGSTAB, and arrow matvec.

All four headers currently embed their own coarse inner/outer iteration logic.

### `MOD_DebugChecks.f90` / `Intermediates`

- reaction-rate calculation helper.
- diagnostic file writer.

The module is compiled but no active caller was found.

## 19. Items already confirmed complete

- The remaining hard-coded Quad convergence `write(2,...)` was already migrated to
  `files%convergence`.
- `Results_Writing` now uses its `n_tot` and `univ` dummy arguments.
- Universes are written once from `univ`, rather than writing the full array inside
  an `n_tot` loop.
- `Header.f90` passes the universe array constructor from `LO_Mesh(:)%univ`.

These completed tasks still need output-contract tests before later refactoring.
