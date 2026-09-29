# kernel/utilities/rlx_split.m

## Purpose

Splits a relaxation superoperator into longitudinal (`R1`), transverse (`R2`) and mixed (`Rm`) components.

## Behaviour

- Syntax: `[R1,R2,Rm]=rlx_split(spin_system,R)`.
- The function first runs a consistency check (`grumble`) that requires:
  - `spin_system` to contain the fields `bas` and `bas.formalism`;
  - the formalism to be `sphten-liouv`;
  - `R` to be a numeric square matrix.
- The spherical tensor basis is interpreted with `lin2lm(spin_system.bas.basis)`, returning rank (`L`) and projection (`M`) arrays.
- Single-spin order states are identified as basis rows whose logical sum over spins equals 1 (`sso_mask`).
- Longitudinal single-spin states satisfy `L>0` and `M==0`; transverse single-spin states satisfy `L>0` and `M~=0`.
- `R1` is a zero matrix of the size of `R` with the block `R(long_sso_mask,long_sso_mask)` copied in; `R2` is built analogously from the transverse block; `Rm` is the remainder `R-R1-R2`.

## Inputs and outputs

**Inputs**

- `spin_system` — spin system object carrying the `sphten-liouv` basis information.
- `R` — relaxation superoperator in `sphten-liouv` formalism; must be a numeric square matrix.

**Outputs**

- `R1` — the part of `R` acting on purely longitudinal single-spin states.
- `R2` — the part of `R` acting on purely transverse single-spin states.
- `Rm` — the rest of `R`.

## References

- Source: [kernel/utilities/rlx_split.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rlx_split.m)
- Spin Dynamics Wiki: [rlx_split.m](https://spindynamics.org/wiki/index.php?title=rlx_split.m)
