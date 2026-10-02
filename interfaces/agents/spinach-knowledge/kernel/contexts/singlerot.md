# kernel/contexts/singlerot.m

- Signature: `[answer,sph_grid] = singlerot(spin_system,pulse_sequence,parameters,assumptions)`
- Source: [`kernel/contexts/singlerot.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/contexts/singlerot.m)
- Existing Wiki page: [`singlerot.m`](https://spindynamics.org/wiki/index.php?title=singlerot.m)

## Role and physical contract

`singlerot` sets up one spinning sample/rotor configuration and calls a user-supplied pulse-sequence function once for each orientation in a spherical integration grid. It combines the spin Hamiltonian, relaxation and kinetics generators, and MAS rotor-phase dynamics. The pulse-sequence function defines the actual experiment and the shape/type of each per-orientation answer; this context does not prescribe an answer matrix shape of its own.

For Liouville formalisms (`sphten-liouv` and `zeeman-liouv`), the context expands the spin-space problem by `spc_dim = 2*max_rank+1` Fourier rotor harmonics. For each crystallite orientation it builds the Hamiltonian blocks over those rotor phases and forms a Fokker-Planck generator from the block-diagonal Hamiltonian and rotor derivative. The rotor term is `M = 2*pi*rate*kron(d_dphi,I)`; the returned generator passed to the sequence is cleaned up using `spin_system.tols.liouv_zero`. Relaxation and kinetics operators are lifted across the rotor subspace. The resulting composite problem has dimension `spc_dim*spn_dim`, where `spn_dim` is the Hamiltonian dimension reported by the source.

For Hilbert formalisms (`zeeman-hilb` and `zeeman-wavef`), the context passes the Hamiltonian stack, one matrix per rotor phase, together with the relaxation and kinetics operators. It does not construct the Liouville rotor derivative in this branch. The pulse sequence is called as `pulse_sequence(spin_system,parameters,G,R,K)` in Liouville space or `pulse_sequence(spin_system,parameters,H,R,K)` in Hilbert space; `parameters.spc_dim` and `parameters.spn_dim` are added for the sequence.

## Inputs and conventions

- `spin_system` is the configured Spinach system. This context obtains its Hamiltonian, relaxation, and kinetics generators and uses `spin_system.bas.formalism` to select the branch. It also reads `spin_system.sys.root_dir` to load grid data from `kernel/grids`, `spin_system.comp.isotopes` to validate rotating-frame spins, `spin_system.tols.liouv_zero` for Liouville cleanup, and temporarily changes `spin_system.sys.output` when the sequence is silenced.
- `pulse_sequence` must be a function handle. The source points to the shipped pulse-sequence functions in the `experiments` directory.
- `parameters.rate` is the signed rotor rate in hertz. The source convention is positive for JEOL and negative for Varian and Bruker, reflecting the rotation-direction convention.
- `parameters.axis` is a normalised, three-component row vector giving the rotor-axis direction in Cartesian coordinates. The source checks that it is a real numeric row with three elements and uses its direction to obtain rotor orientation angles.
- `parameters.spins` is a nonempty cell array of active spin labels, for example `{'1H','13C'}`. `parameters.offset` gives transmitter offsets in hertz for those spins, with the same number of entries; if omitted, it defaults to numeric zeros. Supply a numeric array, not a cell array: the executable grumbler rejects nonnumeric offsets despite the stale header wording.
- `parameters.max_rank` is a nonnegative integer setting the highest retained rotor harmonic. The rotor subspace size is `2*max_rank+1`; convergence is to be checked by increasing the rank, with expected sideband count offered by the source as a starting estimate.
- `parameters.grid` names a spherical grid file in the kernel grids directory. It supplies `alphas`, `betas`, `gammas`, and quadrature `weights`. The source specifies two-angle grids for Liouville space and three-angle grids for Hilbert space.
- `parameters.rframes` describes optional rotating-frame transformations. The source example `{{'13C',2},{'14N',3}}` means second order for carbon-13 and third order for nitrogen-14; the respective spins should use laboratory-frame assumptions.
- `parameters.needs` may request `'iso_eq'`, which makes the context place the equilibrium state of the isotropic Hamiltonian in `parameters.rho0`. A user-supplied `rho0` conflicts with that request, and this option is unavailable for the wavefunction formalism.
- `assumptions` is a character-string context passed to Spinach's assumption handling; the pulse-sequence header documents which assumptions it expects.
- Defaults include no decoupling, no additional rotating frames, zero offsets, silenced sequence output, powder summation enabled, and no extra `needs` flags.

## Orientation, initial state, and return values

The grid provides one quadrature orientation per weight. At each orientation, Wigner rotations combine the molecular orientation, the rotor-axis direction in the laboratory frame, and the rotor phase. Rotor phases are generated on a uniform Fourier grid with `2*max_rank+1` points. The Liouville initial state and coil state, when present, are lifted into the rotor subspace: a single-crystal grid starts at the first rotor phase, while powder simulations start with equal rotor-phase weights. A powder coil state is similarly replicated over the rotor phases.

- `answer`: with `parameters.sum_up` enabled (default), the quadrature-weighted sum of the pulse-sequence outputs over grid orientations. With it disabled, a cell array contains the separate output for each orientation. The component dimensions remain those chosen by the pulse-sequence function.
- `sph_grid`: the loaded grid structure, including the angle arrays and weights used in the calculation.

## Source-supported setup pattern

A source-documented parameter pattern is `parameters.spins = {'1H','13C'}`, a signed `parameters.rate` in hertz, a three-component `parameters.axis`, a grid filename appropriate to the formalism, and an experiment function handle as `pulse_sequence`. For an explicit rotating-frame example, the source gives `parameters.rframes = {{'13C',2},{'14N',3}}`. The source does not include a complete concrete pulse-sequence invocation; use a function from `experiments` for that part.

## Optional FFT rotor derivative

With `sys.enable={'polyadic'}`, the Liouville rotor derivative is a product of three sequential polyadic factors: FFT, Fourier multiplication, and inverse FFT. Ordinary `mtimes` composes them in the order `inverse_fft*multiplier*fft`; their tensor product with the spin identity gives the rotor term. The multiplier is the sparse diagonal matrix with entries `1i*[0:max_rank -max_rank:-1].'`, built once on the CPU as an ordinary numeric core. The FFT adjoint is `spc_dim*ifft`, and the inverse-FFT adjoint is `fft/spc_dim`, preserving MATLAB's transform normalisation. The Hamiltonian phase blocks remain explicit; relaxation and kinetics are lifted as before. The callback receives a polyadic generator and must support exponential-action propagation through `step` or `evolution`. This opt-in retains the same phase grid, rotor direction, and state averaging; it does not lower the chosen rotor rank. The default explicit path and Hilbert branch are unchanged.

When GPU execution is requested, the assembled polyadic generator is uploaded on the executing orientation worker before the pulse-sequence callback. Its numeric multiplier and other uploaded factors are then reused; FFT handles operate on their input blocks without captured device data or per-action multiplier construction.
