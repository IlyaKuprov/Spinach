# kernel/utilities/zte.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/zte.m`
- Signature: `projector=zte(spin_system,L,rho,nstates)`
- Total lines: 191

## Purpose

Zero track elimination inspects the first few trajectory steps and removes coordinates whose amplitudes remain below a user-specified tolerance.

## Physical / mathematical content

Screening follows the actual initial state under the supplied Liouvillian. For a horizontal stack, the retained space contains the union of populated coordinates, without mixing the phases of different columns. This is trajectory-level pruning, not a change of spin basis or an orthogonalised Arnoldi/Lanczos construction.

## Numerical / algorithmic content

Each nonzero column is propagated independently with `step`, which scales that column before its reordered Taylor expansion. Only one column trajectory and its amplitude maxima are held at a time, alongside the aggregate row maxima; additional screening storage is O(N) for an N-by-M input, and a sparse stack is never passed wholesale to the dense propagation routine. Exactly zero columns contribute no support.

Sampling uses `1/cheap_norm(L)` (unit time for a zero generator), up to `zte_nsteps` samples including the initial state. Each column stops when the number of coordinates whose sampled maximum exceeds `zte_tol` stops growing. Without `nstates`, coordinates whose maxima are strictly below `zte_tol` are dropped. With `nstates`, the largest row maxima over columns and sampled times determine the retained coordinates.

The existing explicit-disable, occupied-row-density, and small-matrix-1-norm shortcuts take precedence over propagation and over `nstates`; they return scalar `1`, leaving the basis unchanged. The count is validated against the number of rows, not the number of input columns.

## Syntax

```matlab
projector=zte(spin_system,L,rho,nstates)
```

## Parameters / inputs

- `spin_system`: Spinach system in `zeeman-liouv` or `sphten-liouv` formalism, supplying tolerances and algorithm switches.
- `L`: square Liouvillian used for time propagation.
- `rho`: initial state column or horizontal stack of state columns, with the same number of rows as `L`.
- `nstates`: existing optional positive integer, no larger than the state-space dimension. When screening runs, keeps this number of the most populated coordinates irrespective of the amplitude tolerance.

## Outputs

`projector` projects into the reduced space. Use `L_reduced=P'*L*P` and `rho_reduced=P'*rho`. A scalar `1` signals an unchanged basis.

## Header notes

Set `sys.tols.zte_tol` before `create` to change the default tolerance. With tiny interactions or nearly equivalent spins, disable ZTE by adding `'zte'` to `sys.disable`.

Method reference: [Kuprov, JMR (2008), doi:10.1016/j.jmr.2008.08.008](https://doi.org/10.1016/j.jmr.2008.08.008). See also [the function Wiki page](https://spindynamics.org/wiki/index.php?title=zte.m).
