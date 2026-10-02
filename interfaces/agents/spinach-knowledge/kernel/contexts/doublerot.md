# kernel/contexts/doublerot.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/contexts/doublerot.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=doublerot.m)

## Contract

`doublerot(spin_system,pulse_sequence,parameters,assumptions)` is the double-angle-spinning context. The sequence is a function handle, and `assumptions` is passed to `assume` while constructing the Hamiltonian. The context obtains the spin Hamiltonian and its spherical interaction components from the selected Spinach basis and system; it also builds relaxation and kinetics operators and evaluates the sequence over a powder grid.

In Liouville formalisms (`sphten-liouv` and `zeeman-liouv`), it constructs a Fokker–Planck generator for the two independent rotor phases. In Hilbert formalisms (`zeeman-hilb` and `zeeman-wavef`), it instead supplies a stack of spin Hamiltonians indexed by pairs of rotor phases; that route has no rotor-derivative operator. The two rotor phase axes have `2*rank_outer+1` and `2*rank_inner+1` points, so the spatial dimension is their product. The spin dimension is `size(H,1)`; the context passes these dimensions as `parameters.spc_dim` and `parameters.spn_dim`. In the Liouville route the combined generator dimension is their product, and the rotor derivative term uses `2*pi*rate` with rates in hertz.

## Parameters and grids

- `parameters.rate_outer` and `rate_inner`: rotor rates in Hz.
- `parameters.axis_outer` and `axis_inner`: normalised three-component vectors for the rotor axes.
- `parameters.rank_outer` and `rank_inner`: retained harmonic ranks; increasing them increases the corresponding phase-grid size and retained Fourier content.
- `parameters.grid`: spherical orientation-averaging grid. This is distinct from the two rotor-phase grids. The returned `sph_grid` contains Euler angles and quadrature weights.
- `parameters.rframes`, offsets, and sequence-specific fields may further configure the calculation. `parameters.needs` defaults to `{}` and may contain only `'iso_eq'`. Set `parameters.needs={'iso_eq'}` when the context must create isotropic thermal-equilibrium `rho0`; this overwrites a user-supplied `rho0`. Otherwise provide the initial condition required by the chosen pulse sequence.

The powder weights combine orientation-level sequence outputs when `parameters.sum_up` is enabled; with it disabled, the orientation results are returned separately. Liouville calculations accept two-angle spherical grids. The Hilbert-space rotor-stack route requires a three-angle grid when more than one orientation is used. Although the source header cautions about powder state-projector treatment, the implemented Liouville branch supports `parameters.grid='single_crystal'`: it places `rho0` at the first rotor phase, switches off rotor-phase averaging, and runs the sole grid orientation. Do not generalise that branch to the distinct Hilbert rotor-stack route.

## Source-supported example

`examples/nmr_solids/dor_powder_nav_fplanck_time.m` uses the Liouville route for 14N, with outer/inner rates of 1 MHz and 5 MHz, ranks 7 and 4, and the `rep_2ang_100pts_oct` orientation grid. Its header explicitly says those spinning frequencies are intentionally high to shorten this example, not representative experimental settings. The related frequency-domain example is `examples/nmr_solids/dor_powder_nav_fplanck_freq.m`.

## Polyadic Liouville route

With `sys.enable={'polyadic'}`, both phase derivatives are three-factor FFT polyadics. Their Kronecker sum retains outer-phase, inner-phase, and spin ordering, with each rate in Hz multiplied by `2*pi`. Rotor-dependent Hamiltonian blocks and lifted relaxation/kinetics remain polyadic. The callback must support implicit exponential actions such as `step` or `evolution`. GPU factors are uploaded once per executing orientation worker before the callback. The explicit route, rotor ranks, powder projection, and Hilbert rotor stacks retain their previous meaning.
