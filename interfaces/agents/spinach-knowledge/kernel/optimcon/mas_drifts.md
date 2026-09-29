# kernel/optimcon/mas_drifts.m

- Signature: `drifts=mas_drifts(spin_system,parameters)`

## Purpose and returned ensemble

The function constructs second-order rotating-frame drift Hamiltonians for a quadrupolar nucleus under magic-angle spinning. It samples each crystallite orientation on a two-angle powder grid and each requested initial rotor phase, discretising a rotor period into phase ticks. The Hamiltonian is held constant within each tick. The output is an ensemble cell array: grid orientations occupy the outer index, initial rotor phases the inner index, and each member contains `parameters.n_slices` Hamiltonian matrices in pulse order.

## Parameters

| Field | Meaning and source-side constraint |
|---|---|
| `spins` | One-element cell array containing the isotope string for the nucleus (the source example is `{'27Al'}`); it must identify an isotope in `spin_system`. |
| `axis` | Normalised three-element row vector for the spinning axis. |
| `grid` | Name of a two-angle powder grid. Its weights must be uniform: the implementation rejects nonuniform weights because the optimisation ensemble average is unweighted. |
| `n_ticks` | Number of discrete rotor-phase ticks per rotor period. |
| `n_phases` | Number of initial rotor phases; must divide `n_ticks`. These phases are separated by `n_ticks/n_phases` ticks. |
| `n_slices` | Positive integer number of rotor ticks sampled for the pulse; indexing wraps around the rotor stack. |

The caller is expected to set the spin-system laboratory-frame assumptions first (the source explicitly directs use of `assume()`).

## Construction

After consistency checks, the routine obtains the laboratory-frame Hamiltonian and spherical-tensor components from `hamiltonian`, and the selected nucleus's carrier Hamiltonian from `carrier`. For each powder orientation and rotor tick, it composes Wigner rotations for the crystallite orientation, rotor-axis tilt, and instantaneous rotor phase, then adds the rotated tensor contributions to the isotropic Hamiltonian. It symmetrises the result and applies `rotframe(...,2)` for the selected nucleus, i.e. a second-order rotating-frame transformation. The rotor stack is then cyclically shifted to form the requested initial-phase/pulse-slice members. Orientations are processed in a `parfor` loop.

The grid's two angles do not mean the implementation discards a third Euler angle: the source notes that it uses all three Euler angles supplied by the grid, including the azimuth held in the third angle for two-angle grids.

## Source and references

- MATLAB source: [kernel/optimcon/mas_drifts.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/mas_drifts.m)
- Spin Dynamics Wiki: <https://spindynamics.org/wiki/index.php?title=mas_drifts.m>
