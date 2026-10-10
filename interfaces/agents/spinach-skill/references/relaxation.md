# Relaxation, equilibrium, kinetics, and orientational averaging

## Contents

- [Selecting a relaxation theory](#selecting-a-relaxation-theory)
- [Redfield theory](#redfield-theory)
- [User-supplied rates](#user-supplied-rates)
- [Scalar relaxation](#scalar-relaxation)
- [Which terms of R survive](#which-terms-of-r-survive)
- [Thermal equilibrium and temperature](#thermal-equilibrium-and-temperature)
- [Chemical kinetics](#chemical-kinetics)
- [Powder grids](#powder-grids)

Relaxation and chemical kinetics enter the Liouvillian as dissipative
generators:

```matlab
L=H+1i*R+1i*K;
```

`R=relaxation(spin_system,euler_angles)` and `K=kinetics(spin_system)`;
standard contexts supply these terms to their sequence when `K` is a matrix.
A nonzero higher-order network or time-dependent rate produces `K(t,eta)`,
which ordinary sequences such as COSY and NOESY cannot add to matrices.
Use the `step`/`iserstep` generator-handle route, or a custom sequence that
evaluates the handle at the required stage time and state; passing such a
network to a standard linear sequence is not supported. For example, a
chemistry-only step uses `{ @(t,eta)1i*K(t,eta), t, 'RKMK4' }`.
The Euler angles are optional and used only
by theories supporting relaxation anisotropy: `powder` recomputes `R` at
every grid orientation and `crystal` passes `parameters.orientation`, but
`liquid` and `singlerot` call `relaxation` without angles, so anisotropic
rates have no effect there.

For a rigid molecular geometry without measured tumbling data,
`rotcorr(atom_symbols,xyz,solvent,temperature)` provides a surface-ellipsoid
estimate with shipped elemental radii and temperature-dependent pure-water
or chloroform viscosity. `atom_symbols` is an N-by-1 cell array, `xyz` is
N-by-3 in Angstrom, and `solvent` is the character string `'water'` or
`'chloroform'`. Surface sampling refines automatically to 1% successive-grid
self-consistency, not 1% physical accuracy. Water uses a protein-derived
hydration layer; chloroform uses an uncalibrated bare van der Waals envelope.
The scalar output is an isotropic-equivalent rank-2 time, while the tensor
retains anisotropy in the input frame. Do not interpret the scalar as every
anisotropic correlation time or transfer the protein model's accuracy to
small molecules. Read `etc/rotcorr_sources.md` for parameters, temperature
ranges, sources, and applicability, and the `kernel/utilities/rotcorr`
knowledge entry for usage.

In the sphten direct sum, `unit_state` populates every substance unit coordinate,
and `equilibrium` returns independently normalised, unweighted blocks. IME acts
on each block separately. T1/T2 rates use the local descriptor columns and
`chem.parts` spin map; diagonal retention and damping preserve every unit.

## Selecting a relaxation theory

`inter.relaxation` is a cell array of strings and more than one may be
given. The accepted set is exactly `'damp'`, `'t1_t2'`, `'redfield'`,
`'naka-zwan'`, `'lindblad'`, `'nottingham'`, `'weizmann'`, `'SRFK'`,
`'SRSK'`; anything else is refused by `create`. If the field is absent, `R` is zero only when there are also no dissipative
bosonic modes. Mode damping or dephasing can contribute automatically; see
`kernel/relaxation.m` and `kernel/utilities/rlx_modes.m`.

| Theory | Physical model | Required inputs | Typical use |
|---|---|---|---|
| `redfield` | Bloch-Redfield-Wangsness theory: rotational modulation of anisotropic interactions | `inter.tau_c` | Liquid-state NMR from first principles: NOE, cross-correlation, DD/CSA interference |
| `t1_t2` | Phenomenological per-spin R1 and R2, optionally anisotropic | `inter.r1_rates`, `inter.r2_rates` | Line widths from experiment; solids; wherever measured rates beat computed ones |
| `damp` | Non-selective damping of every state at one rate | `inter.damp_rate` | Crude broadening, absorbing boundaries, numerical stabilisation |
| `lindblad` | Lindblad-form dissipator from per-spin R1 and R2 | `inter.lind_r1_rates`, `inter.lind_r2_rates` | Strongly driven problems where Redfield theory is invalid |
| spin-phonon (`rlx_phonon.m`, not yet an `inter.relaxation` theory) | Generalised Lindblad dissipator of Saito and Miyashita: a coupling operator, a phonon spectral density `I0*w^alpha`, and a bath temperature, evaluated in the eigenbasis of the current Hamiltonian; returned either as the dressed coupling operator (`'hilb'`) or as the Liouville-space superoperator (`'liouv'`) | arguments of `rlx_phonon(spin_system,H,X,I0,alpha,T,form)`; `experiments/pulsed_field.m` requests the `'hilb'` form on every stair of a field profile and applies the dissipator as Hilbert-space matrix products | Pulsed-field magnetometry of molecular magnets; see `examples/giant_spin/case_studies` |
| `weizmann` | DNP model: electron and nuclear R1/R2 plus dipolar cross-relaxation | `inter.weiz_r1e`, `weiz_r2e`, `weiz_r1n`, `weiz_r2n`, `weiz_r1d`, `weiz_r2d` | Solid-effect and cross-effect DNP |
| `nottingham` | DNP model on the four-level two-electron manifold | `inter.nott_r1e`, `nott_r2e`, `nott_r1n`, `nott_r2n` | Cross-effect DNP; exactly two electrons |
| `SRFK` | Scalar relaxation of the first kind: stochastic modulation of a coupling | `inter.srfk_tau_c`, `inter.srfk_mdepth` | Exchange fast enough to modulate J; conformational averaging |
| `SRSK` | Scalar relaxation of the second kind, Abragam's expressions | `inter.srsk_sources` | Quadrupolar neighbours (14N, 35Cl, 79Br) broadening their partners |
| `naka-zwan` | Nakajima-Zwanzig evaluation of the same rotational-modulation kernel as `redfield`; the two are mutually exclusive | `inter.tau_c`, `inter.nz_shift`, `inter.nz_onshell` | Relaxation beyond the Redfield evaluation of the memory kernel |

Nottingham requires one substance with exactly two electrons. `relaxation`
raises `Spinach:relaxation:nottinghamSubstance` for every segmented descriptor,
even when every substance contains an electron pair. The `create` restriction
of two electrons overall is unchanged. Do not interpret Nottingham nuclear
rate parameters as a standalone nucleus-only model.

The trajectory-integral utility `ngce` supports a single chemical substance.
Segmented inputs raise `Spinach:ngce:segmentedSubstances`, with or without
regularisation; its scalar unit-state projection must not be applied to a
direct sum of substances.

Two settings become mandatory the moment `inter.relaxation` is present:

```matlab
inter.rlx_keep='kite';       % which terms of R survive
inter.equilibrium='zero';    % where relaxation drives the system
```

Formalism restrictions are enforced. `t1_t2`, `redfield`, `naka-zwan`, and
`SRSK` are `sphten-liouv` only. `lindblad`, `nottingham`, `weizmann`, `SRFK`,
and IME thermalisation all require Liouville space (`sphten-liouv` or
`zeeman-liouv`). Only `damp` works in `zeeman-hilb`.

Terms accumulate in a fixed order: `t1_t2`, `redfield`, `naka-zwan`,
`lindblad`, `weizmann`, `nottingham`, `SRFK`, `SRSK`; then the dynamic frequency shift
policy applies, then `inter.rlx_keep` truncates, then `damp` is added, then
thermalisation. Hence `SRSK` reads the superoperator accumulated so far to
extract source spin T1 and T2 and is useless without a companion theory,
while `damp` survives even `inter.rlx_keep='diagonal'` in supported formalisms.
For damp-only models, use `labframe` retention in cross-formalism calculations
(as in `thermal_equilibrium_4` and `thermal_equilibrium_5`): the pre-damping
superoperator is zero, so full retention preserves the same generator without
requesting the unimplemented Zeeman diagonal policy.

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
damping does not destroy the trace; Hilbert space does not. Liouville-space
sequences called with `zeeman-hilb` inputs go through `sim2liouv`, which
projects the unit state out of the converted relaxation superoperator (the
report line reads `unit state exempted from the projected relaxation
superoperator`), so the converted operator equals the native Liouville
`damp` and the trace is conserved on that route as well.

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
SRSK with `zeeman-liouv` is explicitly refused as not implemented; no
alternative high-spin Lindblad model is substituted.
The additive SRSK contribution is built without thermalisation; the chosen
IME or DiBari-Levitt method is applied once to the accumulated spin relaxation.
The recursive contribution excludes mode dissipation; the original-temperature
bosonic dissipators are appended once, after outer spin thermalisation.

## Which terms of R survive

`inter.rlx_keep` is mandatory whenever relaxation is switched on.

| Value | Kept | Notes |
|---|---|---|
| `'diagonal'` | Self-relaxation only | Cheapest; no NOE, no cross-correlation. Unit state protected in `sphten-liouv`; not implemented in `zeeman-liouv` |
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

`equilibrium` checks each Liouville Hamiltonian block on its own unit state;
a vanishing action raises `Spinach:equilibrium:notLeftProduct` with the substance
number. Exactly zero blocks in segmented spherical-tensor systems are exempt,
including spinful zero-Hamiltonian substances; they retain the local unit state.
Any cross-substance entries in the Hamiltonian assembled from `I` and oriented
`Q` raise `Spinach:equilibrium:crossSubstanceHamiltonian` before propagation.

Both `'IME'` and `'dibari'` require `inter.temperature`. IME needs the unit
state population to be exactly 1; general propagation does not enforce initial
normalisation, so a badly normalised state gives incorrect source amplitudes.
IME requires block-diagonal relaxation: cross-substance entries of `R` raise
`Spinach:thermalize:crossSubstanceRelaxation` even when they preserve unit states.

In segmented `sphten-liouv`, `steady` initialises and pins the unit coordinate
`bas.offsets(n)+1` of every substance to one; supplied guesses must obey the
same unweighted normalisation. Both Newton and squaring methods accept this layout.
Each unit column must drive an active coordinate in its own substance block;
`Spinach:steady:unthermalisedSubstance` names any block that fails this check.
Steady-state GRAPE dressing likewise removes every conserved unit direction from
the adjoint solve and requires a traceless target in each substance.
`magpump` likewise sources each target block through its own unit coordinate,
including targets spanning several substances, and rejects unit-state pumping
in every block.

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

Reaction records currently require `sphten-liouv`. Declare each species and
each directed first-order exchange channel explicitly:

```matlab
inter.chem.parts={[1 2 3 4 5],[6 7 8 9 10]};
inter.chem.concs=[20 4];
inter.chem.reactions={struct('reactants',1,'products',2,...
                            'matching',[1 6;2 7;3 8;4 9;5 10],'rate',4),...
                      struct('reactants',2,'products',1,...
                            'matching',[6 1;7 2;8 3;9 4;10 5],'rate',20)};
```

`inter.chem.parts` contains disjoint numeric vectors of global spin indices;
empty entries are spin-free substances. Matching pairs must have identical
isotopes, but different species need not have identical spin counts or basis
topologies. Unmatched source spins are traced out; unmatched product spins
arrive at identity. Missing source descriptors are reported rather than
invented. Couplings across species are rejected by `create`; relaxation
parameters such as `inter.tau_c` remain per substance.

The example is at exchange equilibrium because forward and reverse event
fluxes are equal. `state(spin_system,'Lz','1H')` applies initial concentrations
without a chemistry qualifier; `coil_state(spin_system,'Lz','1H','exact')`
is unweighted detection. `equilibrate(K,c0)` may still calculate stationary
concentrations of the separate classical network `dc/dt=K*c`; the classical
matrix is not a chemistry input field. First-order reaction rates have units
of inverse seconds. Higher-order rates multiply the other reactant
concentrations, so their units depend on the chosen concentration unit.
`kinetics` returns a sparse matrix for constant first-order records, or
`K(t,eta)` for mass action and time-dependent rates. It reads instantaneous
concentrations from unit coordinates, including initially empty products.

For a compact exchange-plus-quadrupolar-relaxation demonstration, see
`examples/kinetics/uf6_collisions.m`. It sweeps the forward collision
rate while holding the reverse lifetime and rotational correlation time
fixed, recomputes stationary populations, and uses chemical-population-
weighted excitation. The full laboratory-frame `slowpass` calculation
leaves ZTE off by default and plots individually peak-normalised spectra.
Its representative F/U pair and phenomenological distorted state are not
a full UF6 molecule or calibrated liquid-state prediction. Keep collision
frequency, distorted-state lifetime, and rotational correlation time
separate when adapting this model; scaling both exchange rates preserves
populations but is a different physical scan.

Spin replacement between molecules uses matching rather than a flux matrix.
For a two-proton molecule exchanging its first proton with a one-proton pool:

```matlab
inter.chem.parts={1:2,3}; inter.chem.concs=[1 1];
inter.chem.reactions={struct('reactants',[1 2],'products',[1 2],...
                            'matching',[1 3;2 2;3 1],'rate',2,...
                            'closure','additive')};
```

Both populations are invariant because the reactant and product stoichiometries
coincide. The additive arrival keeps each reactant's internal orders but not
cross-reactant polarisation products: the retained molecular spin survives,
while correlations involving the departing spin are lost. Intramolecular
permutation records instead transport the mapped multi-spin orders. These are
physical matching models, not an automatic translation of every retired flux
matrix. `test_cwdm_flux` supplies an analytic replacement check. Freezing an
additive generator is justified only when its concentrations and rates remain
constant, not for general mass-action or product-closure networks.

For two one-spin substances, the additive event with `parts={1,2}`,
`concs=[2000 500]`, `reactants=[1 2]`, `products=[1 2]`,
`matching=[1 2;2 1]`, and rate 1 reproduces directed rates 500 and 2000
inverse seconds. An order-m rate constant has units concentration^(1-m)/s.

Radical-pair loss uses one first-order reaction record per channel. For a
single species whose first two spins are electrons:

```matlab
inter.chem.parts={1:3}; inter.chem.concs=1;
inter.chem.reactions={struct('reactants',1,'products',[],...
                            'matching',zeros(0,2),'rate',1e6,...
                            'selector',{{'singlet',[1 2]}}),...
                      struct('reactants',1,'products',[],...
                            'matching',zeros(0,2),'rate',1e5,...
                            'selector',{{'triplet',[1 2]}})};
```

The third spin may be a nucleus. These selectors give Haberkorn loss via the
projector anticommutators. Select `jones-hore-singlet` and
`jones-hore-triplet` for Jones–Hore loss; each channel then subtracts the
complementary projected density from the full density. For nonselective
exponential loss, omit selectors and use one empty-product record with the
sum of the two rates. This choice changes the physics, not only numerics.
To track products, declare their parts and matching: selective arrival uses
the projected source before tracing unmatched spins, as tested by
`test_cwdm_selectors`. Use an explicit scalar `inter.nz_shift` when required;
a general reaction network does not determine a unique scalar lifetime.

For the legacy exponential scalar use the summed channel rates; the legacy
Haberkorn/Jones-Hore scalar approximation uses half their sum. Retired chemistry
fields are rejected rather than accepted alongside reaction records.

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
