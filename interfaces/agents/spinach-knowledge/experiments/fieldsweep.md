# experiments/fieldsweep.m

- Signature: `[spec,parameters]=fieldsweep(spin_system,parameters)`

## Purpose

Compute field-swept powder EPR spectra using an expensive eigenfields algorithm, an explicit spherical grid, and a hard-coded Lorentzian line shape.

## Physical / mathematical content

- The specified spin couples to the microwave field. `parameters.mw_freq` sets the microwave frequency; `parameters.fwhm` sets the Lorentzian full width at half maximum.
- The magnetic field in `sys.magnet` must be set to 1 Tesla, irrespective of the sweep window.

## Numerical / algorithmic content

- The initial spherical grid supplies orientations whose convex hull defines triangles for powder integration. A non-symmetric grid helps avoid transition degeneracies.
- Eigenfields are calculated asynchronously at grid vertices. Triangle contributions are evaluated with a recursive Voitlander integrator and summed to form the spectrum.
- The peak-position tolerance is set to one quarter of a field-axis grid interval.

## Parameters / inputs

- `parameters.grid` — initial spherical grid, ideally non-symmetric to avoid transition degeneracies; `'rep_2ang_100pts_sph'` is a good starting point.
- `parameters.spins` — one-element cell array identifying the spin coupled to the microwave field, e.g. `{'E'}`.
- `parameters.mw_freq` — microwave frequency, Hz.
- `parameters.fwhm` — Lorentzian line FWHM, Tesla.
- `parameters.window` — field sweep window in Tesla, `[Bmin Bmax]`.
- `parameters.npoints` — number of points in the sweep.
- `parameters.tm_tol` — relative transition-moment tolerance; 0.01 is a good starting point.
- `parameters.rspt_order` — perturbation-theory order for eigenfields calculation; 2 is a good starting point, or specify `Inf` for exact diagonalisation.
- `parameters.int_tol` — powder-integration tolerance, balancing speed and integration accuracy.

## Outputs

- `spec` — field-swept EPR spectrum.
- `parameters` — updated parameters structure containing the magnetic field axis for plotting in `parameters.b_axis`.

## Implementation structure

- The experiment checks its inputs, constructs the Hamiltonians and microwave operator, and loads the specified spherical grid.
- It computes orientation-dependent eigenfields at grid vertices, integrates over convex-hull triangles, and sums their contributions.
- Call this experiment directly, without a context.

ilya.kuprov@weizmann.ac.il

<https://spindynamics.org/wiki/index.php?title=fieldsweep.m>