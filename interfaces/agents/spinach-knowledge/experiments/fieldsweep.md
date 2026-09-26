# experiments/fieldsweep.m

- Signature: `[spec,parameters]=fieldsweep(spin_system,parameters)`

## Purpose

Field swept powder EPR spectra. A rough implementation with ex- pensive eigenfields algorithm, an explicit spherical grid, and a hard-coded Lorentzian line shape. Syntax: [b_axis,spec]=fieldsweep(spin_system,parameters)

## Physical / mathematical content

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Parameters / inputs

- parameters.grid -initial spherical grid, ideally
- a non-symmetric one to avoid
- transition degeneracies, a good
- start is
- 'rep_2ang_100pts_sph'
- parameters.spins -a one-element cell array speci-
- fying the spin that is coupled
- couple to the microwave field,
- e.g. {'E'}
- parameters.mw_freq -microwave frequency, Hz
- parameters.fwhm -Lorentzian line FWHM, Tesla
- parameters.window -field sweep window in Tesla,
- as a vector [Bmin Bmax]
- parameters.npoints -number of points in the sweep
- parameters.tm_tol -relative transition moment to-
- lerance, 0.01 is a good start
- parameters.rspt_order -perturbation theory order for
- eigenfields calculation, 2 is
- a good start; specify Inf for
- exact diagonalisation
- parameters.int_tol -powder integration tolerance,
- a balance between speed and
- integration accuracy

## Outputs

- b_axis -magnetic field axis for plotting
- spec -field-swept EPR spectrum
- Note: irrespective of the actual sweep extents, the magnetic field
- in sys.magnet should be set to 1 Tesla.
- Note: this experiment should be called directly without a context.

## Implementation structure

- Field swept powder EPR spectra. A rough implementation with ex-
- pensive eigenfields algorithm, an explicit spherical grid, and
- a hard-coded Lorentzian line shape. Syntax:
- [b_axis,spec]=fieldsweep(spin_system,parameters)
- parameters.grid - initial spherical grid, ideally
- a non-symmetric one to avoid
- transition degeneracies, a good
- start is
- 'rep_2ang_100pts_sph'
- parameters.spins - a one-element cell array speci-
- fying the spin that is coupled
- couple to the microwave field,
