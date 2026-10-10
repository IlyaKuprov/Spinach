# kernel/utilities/rlx_t1_t2.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rlx_t1_t2.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rlx_t1_t2.m)

## Purpose

Extended T1/T2 relaxation model that returns the relaxation superoperators separately for the longitudinal and the transverse states.

## Behaviour

- Syntax: `[R1Op,R2Op]=rlx_t1_t2(spin_system,euler_angles)`.
- Calls `grumble` to enforce that `spin_system.bas.formalism` is `'sphten-liouv'`; otherwise it errors with `'this function is only available in sphten-liouv formalism.'`.
- Computes spherical tensor ranks (`L`) and projections (`M`) from each local descriptor `spin_system.bas.basis{n}`.
- Preallocates per-spin `r1_rates` and `r2_rates` column vectors.
- For each spin, the R1 and R2 rate specifications are read from the cell arrays `spin_system.rlx.r1_rates{n}` and `spin_system.rlx.r2_rates{n}`. Each entry may be:
  - A numeric scalar: the rate is assigned directly.
  - A numeric 3x3 tensor: Euler angles must be supplied, otherwise the function errors with `'Euler angles must be specified with anisotropic T1/T2 relaxation theory.'`. The orientation vector is computed as `ort=[0 0 1]*euler2dcm(euler_angles(1),euler_angles(2),euler_angles(3))` (noted in the source as matching `alphas=0` of two-angle grids), and the rate is `ort*current_r*_rate*ort'`.
  - A function handle: Euler angles must be supplied (same error otherwise); the rate is obtained by calling the handle as `current_r*_rate(euler_angles(1),euler_angles(2),euler_angles(3))`.
  - Any other specification triggers an error (`'unknown R1 rate specification.'`).
- Euler angles use the ZYZ active convention in radians and specify the system orientation relative to the input orientation; they have no effect when R1 and R2 rates are scalars.
- After filling the rates, the function verifies that all R1 and R2 rates are real, erroring with `'all R1 and R2 relaxation rates must be real numbers.'` if not.
- Builds diagonal superoperators over the direct-sum Liouville space of dimension `bas.offsets(end)`, mapping local descriptor columns through `chem.parts{n}` and placing rates at the block offsets:
  - Spins in the unit state (`L(n,:)==0`) do not contribute.
  - Spins in longitudinal states (`M(n,:)==0` among contributing spins) contribute their R1 rate.
  - Spins in transverse states (`M(n,:)~=0` among contributing spins) contribute their R2 rate.
  - Each state's total R1 and R2 rates are the sums of the contributing single-spin rates; multi-spin orders relax at the sum of the rates of their constituent single-spin orders.
- Returns `R1Op=-spdiags(r1_diagonal,0,matrix_dim,matrix_dim)` and `R2Op=-spdiags(r2_diagonal,0,matrix_dim,matrix_dim)` as sparse diagonal superoperators with negated (dissipative) diagonals.

## Inputs and outputs

**Inputs**

- `spin_system` — Spinach spin system object; must use the `sphten-liouv` formalism, with relaxation specifications in `spin_system.rlx.r1_rates` and `spin_system.rlx.r2_rates`.
- `euler_angles` — Three Euler angles (ZYZ active convention, radians) specifying system orientation relative to the input orientation; required when R1 and/or R2 rates are given as 3x3 tensors or function handles, and ignored for scalar rates.

**Outputs**

- `R1Op` — Relaxation superoperator containing all longitudinal relaxation terms.
- `R2Op` — Relaxation superoperator containing all transverse relaxation terms.

## References

- Spinach Wiki: [rlx_t1_t2.m](https://spindynamics.org/wiki/index.php?title=rlx_t1_t2.m)
