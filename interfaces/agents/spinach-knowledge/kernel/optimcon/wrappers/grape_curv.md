# kernel/optimcon/wrappers/grape_curv.m

- Source: `kernel/optimcon/wrappers/grape_curv.m`
- Signature: `[traj_data,fidelity,df_du]=grape_curv(waveform_u,u2x,dx_du,spin_system)`
- Total lines: 143

## Purpose

Evaluates a GRAPE objective and its gradient for a waveform expressed in user-defined curvilinear coordinates. The map defines the coefficients of the physical control operators; the optimiser can therefore work in coordinates such as amplitude and phase rather than Cartesian channels.

## Physical / mathematical content

For each time sample, `x=u2x(u)`. The supplied `dx_du(u)` has curvilinear coordinates along rows and Cartesian controls along columns: entry `(j,k)` is `dx(k)/du(j)`. Thus the gradient is `df_du=dx_du(u)*df_dx`, including every Cartesian contribution even if some curvilinear inputs are frozen. The number of curvilinear coordinates need not equal the number of control operators.

## Numerical / algorithmic content

`grape_xy` supplies the Cartesian objective and derivatives, including ensemble effects and rectilinear penalty terms. `control.freeze` must have exactly the shape of `waveform_u`, or be empty. It is withheld from the Cartesian calculation and applied to each curvilinear gradient channel after the pullback. It changes neither the propagated waveform nor the objective values. Combining a non-empty freeze mask with a waveform basis remains unsupported. This wrapper returns no Hessian.

## Syntax

`[traj_data,fidelity]=grape_curv(waveform_u,u2x,dx_du,spin_system)`

`[traj_data,fidelity,df_du]=grape_curv(waveform_u,u2x,dx_du,spin_system)`

## Parameters / inputs

- `waveform_u`: real array with coordinates in rows and time samples in columns.
- `u2x`: function handle mapping one curvilinear coordinate column to a column of coefficients for the control operators.
- `dx_du`: function handle returning the curvilinear-by-Cartesian Jacobian described above for one coordinate column.
- `spin_system`: Spinach problem configured by `optimcon`.

## Outputs

- `traj_data`: system trajectory data for visualisation and progress reports.
- `fidelity`: state-overlap objective followed by the separately weighted penalty values when present.
- `df_du`: gradient in the input coordinate layout; a separate third-dimension channel is returned for each fidelity or penalty term. Frozen entries are zero in every channel.

## Header notes

Penalties are evaluated in the rectilinear representation, then differentiated through the same coordinate map.
