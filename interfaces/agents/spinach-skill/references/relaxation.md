# Relaxation, equilibrium, kinetics, and orientational averaging

Relaxation and chemical kinetics enter the Liouvillian as dissipative
generators:

```matlab
L=H+1i*R+1i*K;
```

`R=relaxation(spin_system,euler_angles)` and `K=kinetics(spin_system)`;
contexts assemble this for you. The Euler angles are optional and used only
by theories supporting relaxation anisotropy: `powder` recomputes `R` at
every grid orientation and `crystal` passes `parameters.orientation`, but
`liquid` and `singlerot` call `relaxation` without angles, so anisotropic
rates have no effect there.

## Selecting a relaxation theory

`inter.relaxation` is a cell array of strings and more than one may be
given. The accepted set is exactly `'damp'`, `'t1_t2'`, `'redfield'`,
`'naka-zwan'`, `'lindblad'`, `'nottingham'`, `'weizmann'`, `'SRFK'`,
`'SRSK'`; anything else is refused by `create`. If the field is absent,
`R` is zero.

| Theory | Physical model | Required inputs | Typical use |
|---|---|---|---|
| `redfield` | Bloch-Redfield-Wangsness theory: rotational modulation of anisotropic interactions | `inter.tau_c` | Liquid-state NMR from first principles: NOE, cross-correlation, DD/CSA interference |
| `t1_t2` | Phenomenological per-spin R1 and R2, optionally anisotropic | `inter.r1_rates`, `inter.r2_rates` | Line widths from experiment; solids; wherever measured rates beat computed ones |
| `damp` | Non-selective damping of every state at one rate | `inter.damp_rate` | Crude broadening, absorbing boundaries, numerical stabilisation |
| `lindblad` | Lindblad-form dissipator from per-spin R1 and R2 | `inter.lind_r1_rates`, `inter.lind_r2_rates` | Strongly driven problems where Redfield theory is invalid |
| `weizmann` | DNP model: electron and nuclear R1/R2 plus dipolar cross-relaxation | `inter.weiz_r1e`, `weiz_r2e`, `weiz_r1n`, `weiz_r2n`, `weiz_r1d`, `weiz_r2d` | Solid-effect and cross-effect DNP |
| `nottingham` | DNP model on the four-level two-electron manifold | `inter.nott_r1e`, `nott_r2e`, `nott_r1n`, `nott_r2n` | Cross-effect DNP; exactly two electrons |
| `SRFK` | Scalar relaxation of the first kind: stochastic modulation of a coupling | `inter.srfk_tau_c`, `inter.srfk_mdepth` | Exchange fast enough to modulate J; conformational averaging |
| `SRSK` | Scalar relaxation of the second kind, Abragam's expressions | `inter.srsk_sources` | Quadrupolar neighbours (14N, 35Cl, 79Br) broadening their partners |
| `naka-zwan` | Nakajima-Zwanzig evaluation of the same rotational-modulation kernel as `redfield`; the two are mutually exclusive | `inter.tau_c`, `inter.nz_shift`, `inter.nz_onshell` | Relaxation beyond the Redfield evaluation of the memory kernel |

Two settings become mandatory the moment `inter.relaxation` is present:

```matlab
inter.rlx_keep='kite';       % which terms of R survive
inter.equilibrium='zero';    % where relaxation drives the system
```

Formalism restrictions are enforced. `t1_t2`, `redfield`, and
`naka-zwan` are `sphten-liouv` only. `lindblad`, `nottingham`, `weizmann`, `SRFK`, `SRSK`,
and IME thermalisation all require Liouville space (`sphten-liouv` or
`zeeman-liouv`). Only `damp` works in `zeeman-hilb`.

Terms accumulate in a fixed order: `t1_t2`, `redfield`, `naka-zwan`,
`lindblad`, `weizmann`, `nottingham`, `SRFK`, `SRSK`; then the dynamic frequency shift
policy applies, then `inter.rlx_keep` truncates, then `damp` is added, then
thermalisation. Hence `SRSK` reads the superoperator accumulated so far to
extract source spin T1 and T2 and is useless without a companion theory,
while `damp` survives even `inter.rlx_keep='diagonal'`.

## Redfield theory

```matlab
inter.relaxation={'redfield'};
inter.tau_c={200e-12};
inter.rlx_keep='kite';
inter.equilibrium='zero';
```

`inter.tau_c` is a cell array with one element per chemical species declared
in `inter.chem.parts`, each a non-negative real vector of one, two, or three
correlation times in seconds: one for isotropic rotational diffusion, two
for axial diffusion (around and perpendicular to the main axis), three for
rhombic diffusion (about the XX, YY, and ZZ directions of the diffusion
tensor). Zero correlation times are refused.

The laboratory frame Hamiltonian is built internally through
`assume(spin_system,'labframe')` and the Redfield expression integrated
numerically, so the relaxation you get is whatever the anisotropic part of
your Hamiltonian supports: dipolar terms need `inter.coordinates`, CSA needs
`inter.zeeman.matrix` or `eigs`, quadrupolar relaxation needs the
quadrupolar tensors. Isotropic shifts and scalar couplings give none.

Validity is checked rather than assumed: if `1/max(abs(diag(R)))` falls
below any correlation time, the run stops with `T1,2>>tau_c validity
condition violation in Redfield theory`. That is a physics error, not a
numerical one; move to `lindblad`, to `gridfree` with the stochastic
Liouville equation, or to measured rates under `t1_t2`. Evaluation is
parallel by default; `'asyredf'` in `sys.disable` forces the serial path.

## User-supplied rates

Prefer measured rates whenever they exist: Redfield theory reproduces only
the mechanisms your Hamiltonian contains, so paramagnetic impurities,
spin-rotation, exchange, and unmodelled internal motion are invisible to it
and the line widths come out too narrow.

```matlab
inter.relaxation={'t1_t2'};
inter.r1_rates={1.0 1.0 0.5};
inter.r2_rates={5.0 5.0 2.0};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
```

`inter.r1_rates` and `inter.r2_rates` are cell arrays with one element per
spin: a real scalar rate in hertz, a real symmetric 3x3 tensor in hertz, or
a function handle of three Euler angles that must be 2*pi-periodic in all
arguments. Tensor and handle specifications need the orientation, so they
only do anything in a context that passes Euler angles into `relaxation`;
the projection is `ort=[0 0 1]*euler2dcm(alpha,beta,gamma)`, matching the
`alphas=0` convention of the two-angle grids. Multi-spin orders relax at the
sum of the rates of their constituent single-spin orders.

`inter.lind_r1_rates` and `inter.lind_r2_rates` are plain vectors, one
non-negative real entry per spin in hertz. Prefer Lindblad over `t1_t2` when
the system is driven hard enough that the semigroup structure matters: the
dissipator is completely positive by construction. `R2>=R1/2` is enforced on
thermodynamic grounds here and on every Weizmann and Nottingham rate pair.
The Weizmann matrices `inter.weiz_r1d` and `inter.weiz_r2d` are
nspins-by-nspins and non-negative, carrying longitudinal flip-flop and
transverse ZZ cross-relaxation respectively.

`inter.damp_rate` is a real scalar or a 3x3 matrix with non-negative
eigenvalues, in hertz. Given Euler angles it is projected onto the
orientation; without them `mean(diag(damp_rate))` is used and a warning
notes the discarded anisotropy. Liouville space spares the unit state, so
damping does not destroy the trace; Hilbert space does not.

## Scalar relaxation

Scalar relaxation of the first kind treats a coupling whose magnitude is
stochastically modulated:

```matlab
inter.relaxation={'redfield','SRFK'};
inter.srfk_tau_c={[1.0 5e-9]};
inter.srfk_mdepth=cell(numel(sys.isotopes));
inter.srfk_mdepth{1,4}=20;    % RMS modulation depth, Hz
```

`inter.srfk_tau_c` is a cell array of two-element `[weight tau_c]` vectors
giving the exponential components of the correlation function.
`inter.srfk_mdepth` is an nspins-by-nspins cell array of empty matrices or
non-negative scalars in hertz, the root mean square modulation depth of the
corresponding coupling, with zero diagonal. The coupling Hamiltonian is
rebuilt with the couplings replaced by their modulation depths and
integrated against the background Hamiltonian; the Redfield validity
condition is retested against each `srfk_tau_c{n}(2)`.

Scalar relaxation of the second kind covers a fast relaxing quadrupolar
neighbour and needs only the source list, `inter.srsk_sources=[4 14]`. The
T1 and T2 of the sources are read off the superoperator built by the
preceding theories, and Abragam's expressions are applied to every
heteronuclear partner using the isotropic part of the coupling tensor
between them. A source that relaxes too slowly for the treatment to hold
stops the run with `SRSK theory is not applicable: source spin relaxation is
too slow`. Contributions are additive and reported in hertz.

## Which terms of R survive

`inter.rlx_keep` is mandatory whenever relaxation is switched on.

| Value | Kept | Notes |
|---|---|---|
| `'diagonal'` | Self-relaxation only | Cheapest; no NOE, no cross-correlation. Unit state protected |
| `'kite'` | Self-relaxation plus longitudinal cross-relaxation | The NOE-capable minimum and usual liquid-state choice. `sphten-liouv` only |
| `'secular'` | All terms connecting states of equal Zeeman frequency | Secular with respect to the Zeeman Hamiltonian. `sphten-liouv` only |
| `'labframe'` | Everything | Only correct for laboratory frame simulations |

`inter.rlx_dfs` decides the fate of dynamic frequency shifts: `'keep'` or
`'ignore'`, defaulting to `'ignore'`, which takes `R=real(R)`. Keeping them
matters where the imaginary part of R shifts lines measurably, mostly a
paramagnetic and quadrupolar concern.

## Thermal equilibrium and temperature

`inter.temperature` is in kelvin. If absent, `create` assumes 298 K and
prints a warning. Negative and complex values pass input validation, which
is how inverted spin temperatures are specified.

| `inter.equilibrium` | Meaning |
|---|---|
| `'zero'` | Relaxation drives the state vector to zero. Correct wherever the equilibrium magnetisation is subtracted out |
| `'IME'` | Inhomogeneous master equation: `equilibrium.m` supplies the lab frame equilibrium state and R is corrected to drive the system there |
| `'dibari'` | DiBari-Levitt: R is multiplied by the imaginary-time propagator of the lab frame Hamiltonian left side product superoperator |

Both `'IME'` and `'dibari'` require `inter.temperature`. IME needs the unit
state population to be exactly 1, which Spinach cannot check, so a badly
normalised initial condition gives a silently wrong steady state.
DiBari-Levitt is more expensive but better behaved in exotic regimes; it
demands a positive real temperature and refuses the high-temperature
approximation (`inter.temperature=0`), and `equilibrium.m` refuses absolute
zero outright since ground states are commonly degenerate. Where Euler
angles are available both methods recompute the equilibrium state at the
current orientation.

Pulse sequences ask for a thermal initial state through `parameters.needs`,
which places the result in `parameters.rho0`: `'rho_eq'` in `liquid`,
`'iso_eq'` in `powder`, `singlerot`, `doublerot`, `floquet`, and `gridfree`,
`'aniso_eq'` in `powder` and `crystal`.

Run `[r1,r2,t1,t2,R]=relaxan(spin_system)` before trusting any of this: it
prints longitudinal and transverse rates and times for every spin, dynamic
frequency shifts dropped. Compare against measurement first.

## Chemical kinetics

Chemical processes are only available in `sphten-liouv`. Declare the species
first:

```matlab
inter.chem.parts={[1 2 3 4 5],[6 7 8 9 10]};
inter.chem.rates=[-4 +20; +4 -20];
inter.chem.concs=[20 4];
```

`inter.chem.parts` is a cell array of disjoint vectors of spin indices, one
per chemical species. When reactions are declared, all parts must have the
same number of spins in matching isotope order, and the basis subspaces on
either side of the reaction arrow must have identical topology; otherwise
`kinetics` stops with `spin systems on either side of the reaction arrow
have different topologies or basis sets`. Couplings across species are
refused by `create`, and `inter.tau_c` must carry one element per part.

`inter.chem.rates` is a real square matrix of side equal to the number of
parts, in hertz, column-stochastic: element `(i,j)` is the rate of
conversion of species `j` into species `i`, the diagonal carries the total
outflow as a negative number, and column sums must vanish, enforced as
conservation of matter. `inter.chem.concs` holds non-negative initial
concentrations, one per part, required whenever rates are given, and is
applied by `state(spin_system,'Lz','1H','chem')`. In the example above the
species are at exchange equilibrium, their concentrations in inverse ratio
to the forward and backward rates. `equilibrate(K,c0)` returns the
equilibrium concentrations of the linear network `dc/dt=K*c`.

Magnetisation transport between individual spins, as opposed to conversion
between whole species, uses a different pair of fields:

```matlab
inter.chem.flux_rate=zeros(2);
inter.chem.flux_rate(1,2)=5e2;    % from spin 1 to spin 2
inter.chem.flux_rate(2,1)=2e3;    % from spin 2 to spin 1
inter.chem.flux_type='intermolecular';
```

The indexing is the transpose of `inter.chem.rates`: `flux_rate(i,j)` is the
rate at which magnetisation moves from spin `i` to spin `j`, in hertz. The
matrix is nspins by nspins and `flux_type` must accompany it.
`'intramolecular'` transports multi-spin orders along with the single-spin
orders and keeps the correlations; `'intermolecular'` damps them instead,
correct when the transported spin lands in a different molecule.

Radical pair recombination is a triple, all three fields required together:

```matlab
inter.chem.rp_theory='haberkorn';       % or 'jones-hore', 'exponential'
inter.chem.rp_electrons=[1 2];
inter.chem.rp_rates=[1e6 1e5];          % [singlet triplet], Hz
```

`'exponential'` removes population uniformly at the sum of the two rates
with no spin selectivity; `'haberkorn'` applies the singlet and triplet
projectors symmetrically from both sides; `'jones-hore'` adds the
cross-terms coupling the two channels. The choice changes the predicted
magnetic field effect, so it is a physical decision, not a numerical one.

## Powder grids

Every context that averages over orientations reads `parameters.grid`, a
string naming a `.mat` file in `kernel/grids`, without the extension. Each
holds `alphas`, `betas`, `gammas`, and `weights`, the weights normalised to
unit sum. Two-angle grids store `alphas` as zeros and sample `betas` and
`gammas`; single-angle grids sample `betas` only.

| Family | Naming | Coverage |
|---|---|---|
| Repulsion, two-angle | `rep_2ang_<N>pts_<sph\|hem\|oct>`, N = 100 to 25600 | Bak-Nielsen repulsion, the workhorse; full sphere, hemisphere, octant |
| Repulsion, single-angle | `rep_1ang_<N>pts`, N = 100 to 6400 | Beta only |
| Repulsion, three-angle | `rep_3ang_<N>pts`, N = 100 to 12800 | All three Euler angles |
| Lebedev, two-angle | `leb_2ang_rank_<L>`, L = 5, 11, 17, ... 131 | Spherical harmonics exact to rank L |
| Lebedev, single-angle | `leb_1ang_rank_<L>`, L = 3, 7, 15, ... 8191 | Beta only |
| Lebedev, three-angle | `leb_3ang_rank_<L>`, L = 5, 11, ... 131 | Wigner functions to rank L |
| Icosahedral | `icos_2ang_<N>pts`, N = 12, 42, 162, ... 163842 | Icosahedron subdivisions |
| Single orientation | `single_crystal` | One point at zero Euler angles |

The point count in the repulsion file names refers to the parent full sphere
grid; the reduced files hold fewer points, so `rep_2ang_400pts_sph` has 400,
`rep_2ang_400pts_hem` has 199, and `rep_2ang_400pts_oct` has 49. Hemisphere
and octant grids are correct only when the orientation dependence of the
observable is symmetric under the corresponding reflections; antisymmetric
interaction components or tensors that do not share principal axes need the
full sphere.

For a quick check `rep_2ang_100pts_sph` shows whether a lineshape is roughly
right; `rep_2ang_800pts_sph` is the commonest choice across the examples and
a sensible default; `rep_2ang_6400pts_sph` and above are for converged
published lineshapes. Lebedev grids are preferable when the spherical rank
of the problem is known, since the file name states the rank integrated
exactly.

Convergence is tested by rerunning on the next grid up and confirming that
the result stops moving. `parameters.sum_up=false` in `powder` returns the
per-orientation outputs instead of the average, exposing grid artefacts
directly, and the second output of `powder`, `singlerot`, `doublerot`, and
`floquet` returns the grid used.
`grid_test(alphas,betas,gammas,weights,ranks,sfun)` plots the residual norm
of spherical function integration against rank, with `sfun` set to `'Y_l0'`,
`'Y_lm'`, or `'D_lmn'` for single-, two-, and three-angle grids. Custom
grids come from `repulsion`, `shrewd`, `grid_fibon`, `grid_igloo`, and
`grid_kron`.

Spinning contexts add a second convergence axis, `parameters.max_rank`,
which must be converged independently of the grid. `gridfree` needs no grid
at all, and takes its rotational correlation times in `parameters.tau_c`
rather than `inter.tau_c`, the latter being read only by Redfield theory.
