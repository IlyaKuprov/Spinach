# Starting points by physical problem

Paths are relative to the Spinach repository root. Every simulation follows the
seven-part shape given in `SKILL.md`; what changes between problem classes is
the context, the assumptions string, the basis, the pulse sequence, and the
processing chain. Skeletons use placeholders of the form `<...>` with the input
unit stated. Fill them from literature, experiment, or a quantum chemistry
calculation, never from plausibility.

## Index of `examples/`

| Directory | Physics |
|---|---|
| `benchmarks` | GPU, parallel and polyadic performance probes |
| `dnp_liq` | Overhauser DNP in liquids, CW microwaves |
| `dnp_mas` | MAS-DNP solid and cross effect |
| `dnp_sol` | Static solid DNP: field/frequency scans, CP, pulsed schemes |
| `esr_liq_pulsed` | Tumbling EPR: FFT pulse-acquire, ENDOR, radical relaxation |
| `esr_sol_pulsed` | Pulse EPR in solids: ESEEM, HYSCORE, DEER, RIDME, SIFTER, ENDOR |
| `esr_sol_swept` | CW and field-swept EPR powder spectra |
| `extremes` | Very large or highly symmetric systems; stress tests |
| `fitting` | Least-squares fitting of experimental spectra |
| `fundamentals` | State-space restriction, symmetry, strong coupling, PFG, spin locking |
| `giant_spin` | Lanthanide and high-spin Hamiltonians, Stevens operators, magnetisation |
| `imaging` | MRI: echoes, phase encoding, EPI, PRESS, spiral, diffusion weighting |
| `karplus_curves` | Karplus coefficients from DFT dihedral scans |
| `kinetics` | Chemical exchange, flux, hyperpolarisation relay, MAS exchange |
| `liquid_crystals` | RDCs from a user-supplied order matrix |
| `microfluidics` | Flow, diffusion, reaction and NMR on a finite-element mesh |
| `nmr_diffusion` | Diffusion, flow and slosh solvers |
| `nmr_liquids` | Pulse-acquire, COSY, HSQC/HMQC/HMBC, NOESY/ROESY, TOCSY, DEPT/INEPT |
| `nmr_metabol` | Metabolite proton spectra from the GISSMO database |
| `nmr_nucleic` | RNA HSQC and NOESY from PDB plus a shift table |
| `nmr_overtone` | Nitrogen-14 overtone NMR (MAS, DOR, CPMAS, DANTE) |
| `nmr_paramag` | Pseudocontact shift, Curie relaxation, DFT spin density |
| `nmr_proteins` | Triple resonance and protein HSQC/NOESY from PDB plus BMRB |
| `nmr_solids` | Static powder patterns, MAS in three formalisms, CP, MQMAS, DOR, REDOR |
| `nmr_spen` | Spatially encoded and ultrafast NMR: UF-COSY, UF-DOSY, PSYCHE |
| `nmr_stochastic` | Primas-style stochastic-excitation NMR; GPU-heavy |
| `nmr_zerofield` | ZULF NMR: zero, Earth, small field, field drop |
| `nqr` | Pure NQR and NQR nutation, powder |
| `optimal_control` | GRAPE, LBFGS and Newton pulse optimisation, state transfer |
| `parahydrogen` | PHIP: PASADENA, ALTADENA, SABRE, ortho-deuterium |
| `quantum_tech` | Spin-cavity QED: Jaynes-Cummings, Purcell, transmons, STIRAP |
| `relaxation_theory` | Redfield superoperators, mechanisms, cross-correlations, TROSY, SLE |
| `shaped_pulses` | Gaussian, chirp, Q5, SLR, VG, fixed-point pulse propagation |
| `singlet_states` | Long-lived singlets: M2S/S2M, decoherence, singlet imaging |
| `spin_chemistry` | Radical-pair singlet yields, yield anisotropy, CIDNP |
| `standard_systems` | Data only: shared Gaussian, ORCA, CASTEP, SpinXML inputs |
| `visualisation` | 3D rendering of shielding, EFG and hyperfine tensors |

## Liquid-state NMR, one dimension

`nmr_liquids/pa_strychnine.m` is the reference for a real molecule: spin system
from the built-in `strychnine()` database, `IK-2` basis with scalar-coupling
connectivity, line widths from a Redfield superoperator rather than apodisation
alone. `fundamentals/strong_coupling.m` is the same experiment on two
hand-entered spins. To move to a new molecule, replace the `[sys,inter]` block
and re-place `parameters.offset` and `parameters.sweep`; everything after
`create` usually survives. `nmr_metabol/molecule_a.m` imports instead with
`gissmo2spinach('<file>.xml',1)`.

## Liquid-state NMR, homonuclear 2D

`nmr_liquids/noesy_strychnine.m`, `cosy90_strychnine.m`, `tocsy_sucrose.m`,
`roesy_strychnine.m`. NOESY starts from thermal equilibrium and so carries
`parameters.needs={'rho_eq'}` with `inter.equilibrium='IME'` and a temperature;
COSY does not. NOESY also needs `bas.space_level=3` rather than `1`, because
cross peaks come from through-space correlations that scalar-coupling
connectivity does not reach. States quadrature reconstruction:

```matlab
fid=liquid(spin_system,@noesy,parameters,'nmr');
fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}});
fid.sin=apodisation(spin_system,fid.sin,{{'sqcos'},{'sqcos'}});
f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1));
f1_sin=real(fftshift(fft(fid.sin,parameters.zerofill(2),1),1));
spectrum=fftshift(fft(f1_cos-1i*f1_sin,parameters.zerofill(1),2),2);
plot_2d(spin_system,-real(spectrum),parameters,20,[0.01 0.1 0.01 0.1],2,256,6,'both');
```

The leading minus is the NOESY sign convention; HSQC-type spectra omit it.

## Liquid-state NMR, heteronuclear 2D

`nmr_liquids/hsqc_strychnine.m` differs in two ways: echo/antiecho rather than
States quadrature, and a split into natural-abundance isotopomers.

```matlab
spin_system=create(sys,inter);
subsystems=dilute(spin_system,'13C');
parfor n=1:numel(subsystems)
    subsystem=basis(subsystems{n},bas);          % basis inside the loop
    fid=liquid(subsystem,@hsqc,parameters,'nmr');
    fid.pos=apodisation(spin_system,fid.pos,{{'sqcos'},{'sqcos'}});
    fid.neg=apodisation(spin_system,fid.neg,{{'sqcos'},{'sqcos'}});
    f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1);
    f1_neg=fftshift(fft(fid.neg,parameters.zerofill(2),1),1);
    spectrum=spectrum+fftshift(fft(f1_pos+conj(f1_neg),parameters.zerofill(1),2),2);
end
```

`basis` must be inside the loop because isotopomers differ in spin count.
`parameters.spins` is ordered `{<F1 channel>,<F2 channel>}` and `parameters.J`
is the working one-bond coupling in hertz that sets the INEPT delays. Siblings:
`hmqc_strychnine.m`, `hmbc_sucrose.m`, `hetcor_strychnine.m`,
`inept_strychnine.m`.

## Proteins and nucleic acids

`nmr_proteins/hsqc_ubiquitin_a.m` imports with
`[sys,inter]=protein('<file>.pdb','<file>.bmrb',options)`, where
`options.pdb_mol` is the model index, `options.noshift='delete'` drops spins
with no assigned shift, and `options.select` (`'backbone-hsqc'` and similar) is
what makes a protein tractable by discarding spins the experiment cannot see.
The basis is `IK-1` with `bas.connectivity='scalar_couplings'`;
`sys.tols.inter_cutoff` (Hz) and `sys.tols.prox_cutoff` (Angstrom) prune weak
couplings. Triple-resonance sequences come in pairs, a minimal four-spin
version and a full-protein version: `hnco_simple.m` and `hnco_ubiquitin.m`,
`hnca_simple.m` and `hnca_gb1.m`, plus `hncoca_simple.m` and `hcanh_simple.m`.
Start from the `_simple` file, confirm peak positions, then swap in
`protein()`. Three-dimensional experiments return four quadrature components
(`fid.pos_pos`, `pos_neg`, `neg_pos`, `neg_neg`), transformed along F3, then
F2, then F1, and plotted with `plot_3d`. RNA is the same pattern with
`nuclacid()`, see `nmr_nucleic/rna_hsqc_theo.m`.

## Solid-state NMR, static powder and MAS

`nmr_solids/static_powder_csa.m` is the shortest complete powder simulation:
anisotropic shielding as eigenvalues plus Euler angles, the `powder` context,
a named grid.

```matlab
inter.zeeman.eigs={[<xx> <yy> <zz>]};    % ppm
inter.zeeman.euler={[<a> <b> <g>]};      % radians
sys.disable={'trajlevel'}; bas.projections=+1;
parameters.grid='rep_2ang_6400pts_sph'; parameters.verbose=0;
fid=powder(spin_system,@acquire,parameters,'nmr');
```

`bas.projections=+1` is safe because a single-quantum spectrum needs one total
projection block; `trajlevel` analysis is meaningless per orientation. Grid
choice is a convergence parameter: refine until the pattern stops moving.

The directory is a deliberate matrix, each system appearing as a
`static_powder_*` file and `_fplanck`, `_floquet`, `_gridfree` MAS variants for
CSA, dipolar, quadrupolar, alanine, glycine, sucrose and tryptophan; compare
`mas_powder_csa_fplanck.m` against `mas_powder_csa_gridfree.m` for two
formalisms on identical physics. `singlerot` (Fokker-Planck) is the default and
adds `parameters.rate` (Hz, sign sets the rotation sense),
`parameters.axis=[1 1 1]` for the magic angle, `parameters.max_rank` and a
Lebedev grid; rank and grid are both convergence parameters and must be raised
together. `gridfree` takes the same three but no grid, `doublerot` handles DOR.
Cross polarisation is `cp_powder_static_nh.m` and
`cp_contact_mas_nh_fplanck.m`; recoupling is `redor_curve.m` and
`pdsd_simple.m`.

## Quadrupolar nuclei and NQR

Quadrupolar coupling is a self-coupling on the diagonal of the coupling cell
array, built by `eeqq2nqi` (`nmr_solids/static_powder_nqi_a.m`):

```matlab
inter.coupling.matrix{n,n}=eeqq2nqi(<e2qQ/h, Hz>,<eta>,<spin quantum number>,...
                                    [<a> <b> <g>]);   % Euler angles, radians
```

Half-integer central-transition work is `mqmas_nqi.m`, under MAS use
`mas_powder_nqi_fplanck.m` and siblings, and `rframe_nqi_dante.m` covers
second-order rotating frames and DANTE excitation. Pure NQR sets `sys.magnet=0`
and relies entirely on the quadrupolar tensor (`nqr/pure_nqr_nitrogen.m`);
nitrogen-14 overtone detection is `nmr_overtone/mas_glycine_1.m`.

## EPR, field-swept and CW

`esr_sol_swept/fieldsweep_nitroxide.m` is the canonical CW powder calculation.
It uses no context: `fieldsweep` is its own driver, works in `zeeman-hilb`, and
`sys.magnet` must be `1` because the real field is in `parameters.window`.

```matlab
sys.magnet=1;                                  % must be 1 for a field sweep
inter.zeeman.matrix{1}=<3x3 g-tensor, dimensionless>;
inter.coupling.matrix{1,2}=<3x3 hyperfine tensor, Hz>;
parameters.mw_freq=<microwave frequency, Hz>;
parameters.fwhm=<line width, tesla>;
parameters.window=[<lower> <upper>];           % tesla
parameters.rspt_order=Inf;                     % exact, not perturbative
parameters.rho0=-state(spin_system,'Lz','E');  % high-temperature approximation
[spec,parameters]=fieldsweep(spin_system,parameters);
```

Hyperfine couplings quoted in gauss must pass through `gauss2mhz` and then be
multiplied by `1e6` to convert MHz to Hz; millitesla values can pass through
`mt2hz` directly. Liquid-phase FFT EPR is
`esr_liq_pulsed/pulse_acquire_methyl.m`, using `liquid` with `'esr'`
assumptions and a common line width. Slow tumbling needs the stochastic
Liouville equation: `relaxation_theory/sle_esr_nitroxide_1.m`.

## EPR, pulsed: ESEEM, HYSCORE, ENDOR

`esr_sol_pulsed/eseem_nitroxide_powder.m` runs `powder` with `@eseem` and
`'esr'` assumptions, starting from `Lz` on the electron, detecting `L+`, with
`parameters.pulse_op` an `Ly` operator and `parameters.screen` the
destination-state screen that prunes the propagation. Apodisation is applied to
`mean(fid)-fid`: subtracting the mean removes the unmodulated component so the
envelope modulation is visible, and this step is specific to ESEEM.
Two-dimensional correlation is `hyscore_nitroxide_powder.m` and echoes are
`hpa_nitroxide_powder.m`. ENDOR files
(`endor_davies_nox_powder.m`, `endor_mims_nox_powder.m`) are brute-force
simulations with explicit soft pulses and orientation selection, and are
expensive; start from `endor_davies_nox_crystal.m` and get line positions right
before paying for a grid.

## DEER

`esr_sol_pulsed/hard_3_pulse_deer_no.m` is the minimal two-label calculation.
The distance enters through `inter.coordinates`, not through a coupling:
Spinach computes the dipolar tensor from geometry, including the orientation
dependence of the g-tensors. The assumptions string is `'deer'`.

```matlab
sys.isotopes={'E','E'};
inter.zeeman.eigs{1}=[<gxx> <gyy> <gzz>];  inter.zeeman.euler{1}=[<a> <b> <g>];
inter.zeeman.eigs{2}=[<gxx> <gyy> <gzz>];  inter.zeeman.euler{2}=[<a> <b> <g>];
inter.coordinates={[0.00 0.00 0.00]; [<x> <y> <z>]};   % Angstrom
parameters.rho0=state(spin_system,'Lz','E');
parameters.coil_prob=state(spin_system,{'L-'},{1});
parameters.ex_prob=operator(spin_system,{'Lx'},{1});
parameters.ex_pump=operator(spin_system,{'Lx'},{2});
parameters.stepsize=<time step, s>; parameters.nsteps=<steps>;
parameters.output='brief'; parameters.grid='rep_2ang_3200pts_sph';
deer=powder(spin_system,@deer_3p_hard_deer,parameters,'deer');
```

The observable is `imag(deer.deer_trace)`, and the basis is `zeeman-hilb` with
no approximation. Variants: `hard_3_pulse_deer_gd_1.m` for other label
chemistries, `hard_3_pulse_deer_exchange.m` when exchange coupling matters,
`soft_4_pulse_deer_2e.m` for shaped pulses and orientation selection,
`ridme_cu_nitroxide.m` and `sifter_nitroxide_powder.m` for the
relaxation-induced and single-frequency alternatives.

## DNP

`dnp_sol/solid_effect_freq_scan_1.m` scans the microwave frequency at one
crystal orientation using the `crystal` context with `@dnp_freq_scan`, an
`IK-0` basis, `parameters.method='lvn-backs'` and
`parameters.needs={'aniso_eq'}`. Irradiation is `parameters.mw_pwr` (rad/s),
`parameters.mw_frq` (rad/s), `parameters.mw_oper` (an `Lx` electron operator)
and `parameters.ez_oper` (`Lz`). Relaxation is the Weizmann model:
`inter.weiz_r1e`, `weiz_r2e`, `weiz_r1n`, `weiz_r2n` in hertz plus the
`weiz_r1d`/`weiz_r2d` dipolar matrices, with `inter.temperature` in kelvin.
Field scans are `solid_effect_field_scan_1.m` and
`cross_effect_field_scan_1.m`; pulsed schemes live in the `top_dnp`,
`xix_dnp` and `novel_dnp` subdirectories.

MAS-DNP has a dedicated non-context driver, `masdnp`, returning the
steady-state enhancement directly (`dnp_mas/solid_effect_mas_powder.m`), with
`*_enlev.m`, `*_dynam.m` and `*_steady.m` companions decomposing the same
physics; these use `inter.equilibrium='dibari'`, the correct treatment for a
strongly polarised electron. Liquid-phase Overhauser DNP is
`dnp_liq/odnp_liquid_1.m`, running `liquid` with `@dnp_time_dep` and `'esr'`
assumptions; transfer comes from Redfield dipolar cross-relaxation, so geometry
and `inter.tau_c` are the physics, and detection is a row of `Lz` states so
each build-up is tracked separately.

## PHIP and SABRE

`parahydrogen/pasadena_propanal.m` is high-field PHIP: the hyperpolarised
initial state is two-spin longitudinal order on the former parahydrogen
protons, `state(spin_system,{'Lz','Lz'},{<spin a>,<spin b>})`, and the sequence
is `@hp_acquire` with a small `parameters.pulse_angle` in radians.
`altadena_propanal.m` is the low-field variant. `sabre_pyridine.m` is the most
instructive file here because it does field cycling explicitly: create at the
polarisation field, start from `singlet(spin_system,<hydride a>,<hydride b>)`,
evolve, disconnect the hydrides with `decouple`, split the Hamiltonian with
`assume(spin_system,'nmr','zeeman')` and `assume(spin_system,'nmr','couplings')`,
then ramp by stepping `Hc+field*Hz` with `step`. Any field-cycling problem can
be built from that pattern.

## Radical pairs, magnetic field effects, CIDNP

`spin_chemistry/singlet_yield_1.m` computes a MARY curve with
`M=liquid(spin_system,@rydmr_exp,parameters,'labframe')`. `sys.magnet` is `1`
for normalisation and the real sweep is `parameters.fields` in tesla, alongside
`parameters.rates` in hertz, `parameters.electrons` and
`parameters.needs={'zeeman_op'}`. The basis carries `bas.projections=0` and
permutation symmetry on equivalent protons via `bas.sym_spins`/`bas.sym_group`,
and `sys.disable={'zte'}` is required because the singlet start state is not
the thermal one. For Haberkorn or Jones-Hore kinetics use `@rydmr` with
`inter.chem.rp_theory`, `rp_electrons` and `rp_rates` (Hz) supplied together;
the two routes must not be mixed. Yield anisotropy uses `powder` with
`parameters.sum_up=0`, returning per-orientation yields plus the grid structure
(`singlet_yield_anisotropy_1.m`). CIDNP is a different mechanism, handled
phenomenologically with `magpump` in `cidnp_pumping_2.m`; the explicit geminate
treatment, with reactant and product Liouville spaces stacked by hand, is
`cidnp_geminate.m`.

## Relaxation studies

`relaxation_theory/t1t2_strychnine.m` is the shortest useful file in the
library: build the system with a Redfield superoperator and call
`relaxan(spin_system)` for a printed rate table. No context, no pulse sequence.
To extract one rate by hand:

```matlab
R=relaxation(spin_system);
Lz=state(spin_system,'Lz','<isotope>'); Lp=state(spin_system,'L+','<isotope>');
R1=-(Lz'*R*Lz)/(Lz'*Lz);  R2=-(Lp'*R*Lp)/(Lp'*Lp);
```

Mechanisms one at a time: `dd_relaxation_1.m`, `dd_csa_xcorr_1.m`,
`quad_relaxation_1.m`, `scalar_relaxation_1.m`. TROSY selection is
`trosy_nh.m`; measurement experiments are `inv_rec_1.m` and
`cpmg_echo_train.m`; beyond Redfield theory, `sle_solid_limit.m`. Correlation
functions from molecular dynamics are in `relaxation_theory/from_md`.

## Chemical kinetics and exchange

`kinetics/exchange_symmetric.m` is the two-site template. Sites are chemical
subsystems, the rate matrix has columns summing to zero, and the initial state
is spread across subsystems with the `'chem'` qualifier.

```matlab
inter.chem.parts={<spin indices of site 1>,<spin indices of site 2>};
inter.chem.rates=[-<k12> <k21>; <k12> -<k21>];   % Hz
inter.chem.concs=[<c1> <c2>];
parameters.rho0=state(spin_system,'L+','<isotope>','chem');
```

Without `'chem'` the starting magnetisation would sit in one site only. Varying
the rates to walk the spectrum through coalescence is the cheapest validation
of an exchange model. `exchange_asymmetric.m` handles unequal populations,
`flux_symmetric.m` irreversible flux, `glucose_exsy_a.m` two-dimensional EXSY,
`relayed_hyperpol.m` polarisation relayed through a reaction, and
`equilibrate(rates,concs)` the equilibrium composition of a rate matrix.

## MRI, imaging, flow and diffusion

`imaging/phase_encoding_2d.m` is the reference; the `imaging` context takes no
assumptions argument. The spatial problem is a box, a voxel count and a
derivative operator; the sample is a set of phantoms Kronecker-multiplied with
spin states and relaxation operators.

```matlab
parameters.ro_grad_amp=<read gradient, T/m>;   parameters.ro_grad_dur=<s>;
parameters.pe_grad_amp=<phase gradient, T/m>;  parameters.pe_grad_dur=<s>;
parameters.t_echo=<echo time, s>; parameters.image_size=[<px> <px>];
parameters.dims=[<x> <y>];        % metres
parameters.npts=[<nx> <ny>];      % voxels
parameters.deriv={'period',3};
[R1,R2]=rlx_t1_t2(spin_system);
parameters.rlx_ph={<R1 phantom>,<R2 phantom>};  parameters.rlx_op={R1,R2};
parameters.rho0_ph={<phantom>}; parameters.rho0_st={state(spin_system,'Lz','1H')};
parameters.coil_ph={<phantom>}; parameters.coil_st={state(spin_system,'L+','1H')};
mri=imaging(spin_system,@phase_enc_2d,parameters);
```

`sys.disable={'pt','krylov'}` is standard here. Variants: `gradient_echo_1d.m`,
`echo_planar_2d.m`, `press_2d_example.m`, `diffusion_weighted_2d.m`.
`nmr_diffusion/diffusion_test_1.m` solves the
diffusion equation with a ghost spin and `sys.magnet=0`, with no spin dynamics
at all; use it to verify a spatial discretisation before adding magnetisation.
`microfluidics/plain_nmr.m` runs the `meshflow` context on an imported
finite-element mesh, coupling flow, diffusion, reaction and NMR. Ultrafast
experiments are `nmr_spen/ufcosy_2spin.m` and `ufdosy_1spin.m`.

## Optimal control

`optimal_control/state_transfer_pro.m` is the full pattern: build the system,
normalise initial and target states, collect control and offset operators,
assemble the drift, hand a separate `control` structure to `optimcon`, optimise
with `fmaxnewton`.

```matlab
rho_init=state(spin_system,{'Lz'},{<spin>}); rho_init=rho_init/norm(full(rho_init),2);
rho_targ=state(spin_system,{'Lz'},{<spin>}); rho_targ=rho_targ/norm(full(rho_targ),2);
H=frqoffset(spin_system,hamiltonian(assume(spin_system,'nmr')),parameters);
control.drifts={{H}};
control.operators={<Lx and Ly operators per channel>};
control.off_ops={<Lz operators per channel>};
control.offsets={<offset ranges per channel, Hz>};
control.rho_init={rho_init};  control.rho_targ={rho_targ};
control.pwr_levels=2*pi*linspace(<lower>,<upper>,<n>);   % rad/s
control.pulse_dt=<slice duration, s>*ones(1,<n slices>);
control.penalties={'NS'};  control.p_weights=<weight>;
control.method='lbfgs';    control.max_iter=<iterations>;
spin_system=optimcon(spin_system,control);
pulse=fmaxnewton(spin_system,@grape_xy,guess);
```

The state normalisation is not cosmetic: the fidelity functional assumes unit
norm. Sweeping `pwr_levels` and `offsets` makes B1 inhomogeneity and
transmitter misplacement part of the optimisation target rather than something
discovered afterwards. Verify by propagating with `shaped_pulse_xy` and taking
`real(rho_targ'*rho)`. The `features_*.m` files demonstrate one concept each
(amplitude constraints, dissipative drift, time-step optimisation, freezing,
keyholes, multiple targets, phase cycling, wave bases); solid-state control is
`static_powder_control.m` and `mas_powder_control.m`. Analytically designed
rather than optimised pulses are propagated in
`shaped_pulses/shaped_pulse_gaussian.m` and its chirp, Q5 and SLR siblings.

## Fitting to experimental data

`fitting/fumarate_global.m` is the template. A driver loads and normalises the
data, sets a guess and `optimset` options, and calls `fminsearch` with a
closure; a local error function unpacks the parameter vector, rebuilds the spin
system from scratch, simulates, resamples onto the experimental axis, and
returns a scalar residual.

```matlab
function err=errfun(axis_expt,spec_expt,params)
sys.output='hush'; sys.disable={'hygiene'};       % mandatory inside the loop
... absorb params, build sys/inter/bas, create, basis, simulate, apodise, FT ...
axis_theo=sweep2ticks(parameters.offset,parameters.sweep,parameters.zerofill);
spec_theo=interp1(axis_theo,spec_theo,axis_expt,'pchip');
err=norm(spec_expt-spec_theo)^2;
```

Four things make this work. `sys.output='hush'` and `sys.disable={'hygiene'}`
stop Spinach flooding the console and repeating consistency checks on every
function evaluation. Linewidth and amplitude are always fitted parameters, not
fixed. The `interp1(...,'pchip')` resampling onto the experimental axis is what
makes a pointwise residual meaningful. Relative weighting between datasets
comes from the pre-fit normalisation constants and nothing else; there is no
regularisation anywhere in the shipped examples. Variants worth copying:
disjoint spectral windows concatenated into one residual and integral-normalised
to their proton counts (`fitting/fluoroalkanes/syn_difluorobutane.m`);
parameters as offsets from literature values, rates fitted in log space
(`fitting/nmr_kinetics/glucose_exsy_b.m`, which needs `'UseParallel'` and an
`if ~isworkernode` guard around the live plot, and ships with `'MaxIter',10` so
it will not converge as supplied).

Two fits need no optimiser. `fitting/nmr_rdc/saupe_example.m` reads coordinates
with `read_pdb_pro`, matches them to measured couplings, and gets the Saupe
order matrix from `S=rdc_fit(isotopes,xyz,rdc_hz)` by pseudoinverse least
squares, with `xyz2rdc` back-calculating for the correlation plot.
`karplus_curves/leu_chi_fit.m` calls
`[A,B,C,sA,sB,sC]=karplus_fit('<log directory>',{[<i> <j> <k> <l>]})` over a DFT
dihedral scan, where the four indices define the dihedral and the coupling is
read between atoms one and four; the returned curve is
`A*cosd(phi).^2+B*cosd(phi)+C`.

## Singlet states

`singlet_states/m2s_example.m` converts magnetisation to singlet order with no
context: the Hamiltonian is built with `hamiltonian(assume(spin_system,'nmr'))`,
`Lx` and `Ly` operators are passed to the kernel function
`m2s(spin_system,H,Cx,Cy,rho0,<J coupling, Hz>,<Zeeman frequency difference,
Hz>)`, and the population is read off by projecting onto
`singlet(spin_system,<spin a>,<spin b>)`. `s2m_example.m` is the reverse.
Lifetimes limited by intramolecular mechanisms are in the `decoherence_*.m`
files; `singlet_imaging_1.m` combines singlet order with the imaging context.

## Zero-field and low-field NMR

`nmr_zerofield/zero_field_methanol.m` sets `sys.magnet=0`, uses `'labframe'`
assumptions, and calls `@zerofield` through the `liquid` context. Chemical
shifts are irrelevant at zero field, so the spectrum is determined entirely by
J-couplings and magnetogyric ratios; permutation symmetry on equivalent groups
is the main cost saving. `parameters.detection='uniaxial'` and
`parameters.flip_angle` describe the readout, and subtracting `mean(fid)`
before apodisation removes the DC component. Related: `zero_field_benzene.m`,
`small_field_acetonitrile.m`, `field_drop_acetonitrile.m`.

## Paramagnetic NMR and partial alignment

`nmr_paramag/simple_pcs_1.m` introduces a point susceptibility centre through
`inter.suscept.chi` (a 3x3 tensor) and `inter.suscept.xyz` (its position in
Angstrom), giving both a pseudocontact shift and Curie relaxation.
`inter.rlx_keep='labframe'` matters here: the Curie mechanism is not secular in
the usual sense and is lost under a tighter truncation. `point_vs_distr.m`
compares the point-dipole approximation against a distributed susceptibility,
`dft_density.m` uses DFT spin densities, and lanthanide-binding case studies
are in `nmr_paramag/calbindin`.

Partial alignment is a separate mechanism. `liquid_crystals/rdc_twospin.m`
supplies `inter.order_matrix` (a dimensionless 3x3 Saupe matrix) alongside
`inter.coordinates`, and requests the RDC machinery with
`parameters.needs={'rdc'}` before calling `@clip_hsqc`. The couplings follow
from geometry and order matrix; do not also add them by hand.

## Giant spin, lanthanides and molecular magnets

Hilbert-space problems with large multiplicities. `E<N>` declares an electron
of multiplicity `N`, so `E4` is S=3/2, `E8` is S=7/2, `E13` is J=6 and `E16` is
J=15/2. No file in `giant_spin` uses a context; each calls an experiment driver
directly, and `sys.magnet` must be `1` whenever the driver sweeps the field.
Zero-field splitting via D and E (`quartet_levels.m`):

```matlab
sys.magnet=1.0;                        % must be 1 for a field scan
inter.zeeman.matrix={<3x3 g-tensor>};
inter.coupling.matrix{1,1}=zfs2mat(<D, Hz>,<E, Hz>,<a>,<b>,<g>);
parameters.fields=[<lower> <upper>];   % tesla
parameters.orientation=[<a> <b> <g>];  % radians
parameters.nstates=<levels to track>;
fieldscan_enlev(spin_system,parameters);
```

Higher-rank crystal fields go in through `inter.giant.coeff`, a per-spin cell
of rank-indexed vectors of length `2k+1` in hertz, with `inter.giant.euler`
giving a rotation per rank. Stevens parameters from a ligand-field calculation
are converted with `icm2hz` then `stev2sph` and rotated with `wigner`
(`dy_lft_single_1.m`). Drivers: `fieldscan_enlev` for Zeeman diagrams,
`fieldscan_magn` for finite-speed sweep magnetisation and hysteresis, `eqmag`
for equilibrium molar magnetisation (`create` and `basis` must be re-run inside
the loop when field or temperature changes), `fieldsweep` for powder EPR
(`lanthanide_powder.m`), `geffect` for the effective g-tensor of a Kramers
doublet. Relaxation work needs `sphten-liouv`: `lanthanide_redfield.m` for
electron T1 and T2 against the axial ZFS parameter, `nuclear_relaxation_1.m`
for nuclear relaxation near a fast-relaxing ion by adiabatic elimination with
`adelim`. Multi-centre clusters use `sys.enable={'sodd'}` for spin-orbit
corrections to the dipolar couplings; exchange couplings enter
`inter.coupling.scalar` in hertz in the NMR convention.

## Tensor visualisation, quantum technology, fundamentals

Interaction tensors are rendered straight from the parser output, without
`g2spinach`: `props=gparse('<file>.log')` (or `oparse`, `c2spinach`), then
`cst_display(props,{'<element>'},<scaling factor>,[],options)` with
`options.style` set to `'ellipsoids'` or `'harmonics'`. `efg_display` and
`hfc_display` take the same arguments, and the empty fourth argument requests
automatic bonding at a 1.6 Angstrom cutoff. The scaling factor is a display
convenience with no physical meaning, and differs by orders of magnitude
between tensor types because their eigenvalues do. Examples: `cst_strychnine.m`, `efg_silicate.m`, `hfc_pyrene.m`.

`quantum_tech/jaynes_cummings_a.m` couples a spin to a cavity mode with a
finite photon-number ladder; the directory also covers Tavis-Cummings, Purcell
decay, vacuum Rabi splitting, transmons, STIRAP and NV centres
(`quantum_tech/diamond_defects`). `fundamentals` is where to look when a method
rather than an experiment is in question: `strong_coupling.m` and
`roof_effect.m` for second-order spectra, `state_spaces_1.m` for
correlation-order dynamics and the bypass route that assembles `L=H+1i*R` by
hand instead of using a context, and the symmetry and state-space restriction
files for basis-truncation behaviour. `extremes/high_symmetry_1.m` needs tens
of cores and over a hundred gigabytes; `nmr_stochastic/snmr_strychnine.m`
requires a GPU.
