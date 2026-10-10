# kernel/contexts/floquet.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/contexts/floquet.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=floquet.m)

## Contract

`floquet(spin_system,pulse_sequence,parameters,assumptions)` is the magic-angle-spinning powder context. The pulse sequence must be a function handle; the assumptions string is passed to `assume` for Hamiltonian construction. The implementation is restricted to the `sphten-liouv` and `zeeman-liouv` formalisms, builds the Hamiltonian and its spherical interaction components in the selected spin basis, then forms a Fourier/Floquet Liouvillian for each powder orientation. Relaxation and kinetics operators are included in the spin-space part.

With cutoff `max_rank`, the spatial Fourier dimension is `2*max_rank+1`, labelled by ranks from `-max_rank` through `max_rank`; `spn_dim=size(H,1)`. The combined Floquet problem dimension is `spc_dim*spn_dim`. The rotor-turning term is constructed as `2*pi*rate` times the harmonic index, with `rate` in Hz. The implementation detects non-empty spherical interaction ranks and warns if the cutoff is below one of them; the header recommends increasing the cutoff until the result converges, roughly in line with the number of spinning sidebands.

## Parameters and orientation grid

- `parameters.rate`: spinning rate in Hz. The source convention is positive for JEOL and negative for Varian and Bruker, reflecting their rotation directions.
- `parameters.axis`: normalised three-component rotor-axis vector.
- `parameters.max_rank`: required Fourier cutoff.
- `parameters.grid`: spherical averaging grid from the kernel grids directory; its Euler angles and weights define the powder orientations, separately from the Floquet harmonic index. The `single_crystal` grid is explicitly rejected by `floquet()`; use `singlerot()` for a single-crystal simulation.
- `parameters.spins` and `parameters.offset`: spin labels and corresponding transmitter offsets in Hz.
- `parameters.sum_up`: return a weighted orientation average when enabled, or a cell array of per-orientation outputs when disabled. The context also adds `spc_dim` and `spn_dim` to the parameter structure passed to the sequence.

Numerical rotating-frame transformations through `parameters.rframes` are not supported by this Floquet implementation. Its supported `parameters.needs` request is `iso_eq`; the source rejects other needs.

When supplied, `parameters.rho0` and `parameters.coil` are projected into the central Floquet harmonic by `kron(P,...)`. The `iso_eq` request replaces a supplied `rho0` with equilibrium from the isotropic lab-frame Hamiltonian; without it, the context does not invent an initial state. Omitted `parameters.decouple` defaults to no decoupling, and `parameters.serial=true` selects serial orientation evaluation.

## Source-supported example

`examples/nmr_solids/mas_powder_gly_floquet.m` calculates a 13C MAS powder spectrum for glycine at 14.1 T. It uses a 2 kHz rate, `max_rank=23`, and the `leb_2ang_rank_23` grid, with 256 time points and a 50 kHz sweep. These are example settings, not universal convergence prescriptions.

## State-dependent chemistry boundary

This context rejects a function handle returned by `kinetics` with `Spinach:floquet:stateDependentKinetics`. Multi-reactant or callback-rate reaction records require a custom pulse sequence using `step`/`iserstep`, rather than static context assembly; see `examples/kinetics/nonlinear/bimolecular_closures.m` and `examples/microfluidics/reacting_flow_nmr.m`. Constant matrix kinetics remain supported.
