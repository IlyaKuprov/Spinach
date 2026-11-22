# Contexts, pulse sequences, and processing

## Contents

- [The context call](#the-context-call)
- [Writing a pulse sequence](#writing-a-pulse-sequence)
- [The `state` and `operator` grammar](#the-state-and-operator-grammar)
- [The experiments library](#the-experiments-library)
- [Acquisition and processing conventions](#acquisition-and-processing-conventions)

The context sits between the spin system and the pulse sequence: it builds
the Hamiltonian, the relaxation and kinetics superoperators, performs
whatever orientational or spatial averaging the physical situation demands,
and calls the pulse sequence once per orientation or once in total. The pulse
sequence is a plain MATLAB function that receives the generators and returns
whatever it computes; the context returns that unchanged, or averaged.

## The context call

`assumptions` is the string passed to `assume` when the Hamiltonian is built;
`imaging` and `meshflow` take no such argument. The grid-averaging contexts
return the grid they used as a second output.

```matlab
answer=liquid(spin_system,pulse_sequence,parameters,assumptions)
answer=crystal(spin_system,pulse_sequence,parameters,assumptions)
answer=gridfree(spin_system,pulse_sequence,parameters,assumptions)
[answer,sph_grid]=powder(spin_system,pulse_sequence,parameters,assumptions)
[answer,sph_grid]=singlerot(spin_system,pulse_sequence,parameters,assumptions)
[answer,sph_grid]=doublerot(spin_system,pulse_sequence,parameters,assumptions)
[answer,sph_grid]=floquet(spin_system,pulse_sequence,parameters,assumptions)
answer=imaging(spin_system,pulse_sequence,parameters)
answer=meshflow(spin_system,pulse_sequence,parameters)
```

Fields common to the spectroscopic contexts - `liquid`, `crystal`, `powder`,
and the spinning ones:

| Field | Meaning |
|---|---|
| `parameters.spins` | cell array of channel isotopes in channel order, e.g. `{'1H','13C'}` |
| `parameters.offset` | transmitter offset in Hz on each spin listed in `parameters.spins`, in the same order |
| `parameters.rframes` | rotating-frame specification, e.g. `{{'13C',2},{'14N',3}}` for second-order in carbon-13 and third-order in nitrogen-14; the assumptions on those spins must then be laboratory frame |
| `parameters.needs` | cell array of strings requesting extra quantities from the context (below) |
| `parameters.rho0`, `parameters.coil` | initial and detection states, usually built with `state` |

`parameters.needs` is how a sequence asks the context for something only the
context can compute. The accepted strings differ by context:

| String | Effect | Available in |
|---|---|---|
| `'zeeman_op'` | lab-frame Zeeman part of the Hamiltonian, placed in `parameters.hzeeman` | `liquid`, `crystal`, `powder` |
| `'rdc'` | processes residual anisotropic couplings from the user order matrix | `liquid` |
| `'rho_eq'` | thermal equilibrium with respect to the isotropic Hamiltonian, placed in `parameters.rho0` | `liquid` |
| `'iso_eq'` | thermal equilibrium of the isotropic Hamiltonian into `parameters.rho0` | `powder`, `singlerot`, `doublerot`, `floquet`, `gridfree` |
| `'aniso_eq'` | equilibrium recomputed from the full anisotropic Hamiltonian at the current orientation | `crystal`, `powder` |

### Context-specific parameters

| Field | Meaning | Contexts |
|---|---|---|
| `.orientation` | row vector of three Euler angles in radians, relative to the input orientation | `crystal` |
| `.grid` | spherical grid file name from `kernel/grids` (`leb_2ang_rank_*`, `rep_2ang_*`, `icos_2ang_*pts`, `single_crystal`) | `powder`, `singlerot`, `doublerot`, `floquet` |
| `.rate` | spinning rate in Hz; positive for JEOL, negative for Varian and Bruker, whose rotation direction differs | `singlerot`, `floquet`, `gridfree` |
| `.axis` | spinning axis as a normalised three-element vector, the direction about which the rotor turns | `singlerot`, `floquet`, `gridfree` |
| `.max_rank` | maximum rotor harmonic (or Wigner D-function) rank retained; increase to convergence, a good first guess is the number of expected sidebands | `singlerot`, `floquet`, `gridfree` |
| `.rate_outer`, `.rate_inner`, `.axis_outer`, `.axis_inner`, `.rank_outer`, `.rank_inner` | the same three quantities for each of the two rotors | `doublerot` |
| `.tau_c` | rotational diffusion correlation times in seconds; a scalar for isotropic diffusion, a 3x3 matrix for anisotropic | `gridfree` |
| `.serial` | true disables automatic parallelisation | `powder`, `singlerot`, `doublerot`, `floquet` |
| `.sum_up` | 1 (default) returns the powder average, 0 returns a cell array of per-orientation answers | `powder`, `singlerot`, `doublerot`, `floquet` |

In `powder`, `parameters.rho0` may be a function handle of the three ZYZ
active Euler angles. In `singlerot`, two-angle grids belong in Liouville
space and three-angle grids in Hilbert space; `singlerot` supports
arbitrary-order rotating-frame corrections and `floquet` does not.
`doublerot` supports the same opt-in FFT differentiation in both Liouville formalisms. It uses independent outer and inner phase axes, keeps their specified ranks, and passes polyadic Hamiltonian/rotor, relaxation, and kinetics generators to the callback. Use implicit exponential actions; the Hilbert Hamiltonian-stack route is unchanged.

For traditional MAS stacks, `rotor_stack` applies its explicit assumptions
to both Hamiltonian construction and numerical `parameters.rframes`,
regardless of prior assumptions on the input object. Numerical frames reject
all spins under `nmr`/`cavity`, electrons under `esr`/`deer`/`deer-zz`/
`spin-phonon`, and spin-half nuclei under `qnmr`. Nuclear frames remain
allowed under the electron-only sets, as do higher-spin nuclei under `qnmr`.
Numerical frames are not implemented for `se_dnp_h+`, `se_dnp_h-`, or
`se_dnp_h0`: these carrier-free solid-effect components are not complete
laboratory Hamiltonians. Their empty-frame rotor stacks remain available.
`gridfree` performs the spherical integration inside the SLE formalism, so
supplying `parameters.grid` is an error; it is Liouville-space only, and its
correlation times go in `parameters.tau_c`, not `inter.tau_c`, which is used
only by the Redfield module.

`imaging` builds the spatial Fokker-Planck generator from `parameters.u`,
`.v`, `.w` (X, Y, Z velocity components at each sample point, m/s),
`parameters.diff` (diffusion coefficient or 3x3 tensor in m^2/s, for the case
where it is the same in every voxel) or the per-voxel components
`parameters.dxx` … `parameters.dzz`, `parameters.dims` (3D box dimensions in
metres), `parameters.npts` (points per dimension), and `parameters.deriv`,
either `{'fourier'}` or `{'period',n}` for n-point central finite differences
with periodic boundary conditions. Three phantoms are mandatory, each a pair
of a spatial intensity cell array and an operator or state cell array:

```matlab
parameters.rlx_ph={Ph1,...,PhN};  parameters.rlx_op={R1,...,RN};
parameters.rho0_ph={Ph1,...,PhN}; parameters.rho0_st={rho1,...,rhoN};
parameters.coil_ph={Ph1,...,PhN}; parameters.coil_st={rho1,...,rhoN};
```

The relaxation pair is `_ph`/`_op`; the state pairs are `_ph`/`_st`, and the
consistency checks enforce those names, which the file header states
correctly. Each `Ph` has the dimension of the voxel grid, and
each `rho` comes from `state`. The direct product order is Z(x)Y(x)X(x)Spin,
corresponding to column-wise vectorisation of a 3D array with dimensions
ordered `[X Y Z]`. `meshflow`, the microfluidics and magnetohydrodynamics
context on an imported mesh, needs five phantom pairs on the same pattern:
`H_ph`/`H_op`, `R_ph`/`R_op`, `K_ph`/`K_op`, `rho0_ph`/`rho0_st`,
`coil_ph`/`coil_st`.

For imaging flow/diffusion with `parameters.deriv={'fourier'}` and polyadics enabled, `hydrodynamics` retains FFT derivative cores rather than dense matrices. These outputs and the inherited `v2fplanck` flow/diffusion generators are action-only: neither `inflate` nor `full` can materialise them. Disable polyadics when a sparse numeric operator is required. Tensor ordering, physical lengths, and products of first derivatives in `v2fplanck` are unchanged, including even-grid Nyquist treatment. The periodic finite-difference option is unchanged.

### What the context hands to the sequence

The three or five positional arguments after `parameters` are not always what
their names suggest, and a sequence written for one context will not always
run under another.

| Context | Call made |
|---|---|
| `liquid`, `crystal`, `gridfree` | `pulse_sequence(spin_system,parameters,H,R,K)` |
| `powder` | `pulse_sequence(spin_system,localpar,H,R,K)`, once per orientation |
| `singlerot` | Liouville space: the third argument is the Fokker-Planck generator; Hilbert space: a stack of Hamiltonians, one per rotor phase |
| `doublerot` | Liouville space: third argument is `L+1i*M`, including both rotor turning terms; `zeeman-hilb`/`zeeman-wavef`: a cell stack of Hamiltonians over the two rotor phases |
| `floquet` | third argument is the Floquet generator |
| `imaging`, `meshflow` | `pulse_sequence(spin_system,parameters,H,R,K,G,F)` |

`G` is the cell array of gradient operators and `F` the diffusion-and-flow
generator. All contexts that build a spatial subspace also set
`parameters.spc_dim` and `parameters.spn_dim`, the spatial and spin subspace
dimensions, before calling the sequence; `crystal` sets `spc_dim` to 1.

For a standalone periodic Fourier Laplacian, call `fourlap_poly(spin_system,npoints,extents)` to honour the polyadic enable switch. The original two-argument `fourlap` remains explicit. The implicit second derivative retains the even-grid Nyquist mode; it is not a substitute for squaring a first derivative whose Nyquist action is zero.

## Writing a pulse sequence

```matlab
function fid=my_sequence(spin_system,parameters,H,R,K)
```

Assemble the Liouvillian yourself, in this form and no other:

```matlab
L=H+1i*R+1i*K;
```

Build pulse and detection operators with `operator` and `state`, propagate
with `evolution`, `step`, or `krylov`, and return the result. Assume nothing:
check with `isfield` and either supply a documented default or fail. Imaging
sequences additionally receive `G` and `F` and assemble the background
generator as `B=H+F+1i*R+1i*K`.

| Function | Signature |
|---|---|
| `evolution` | `answer=evolution(spin_system,L,coil,rho,timestep,nsteps,output,destination)` |
| `step` | `rho=step(spin_system,L,rho,time_step)` |
| `krylov` | `answer=krylov(spin_system,L,coil,rho,time_step,nsteps,output)` |
| `propagator` | `P=propagator(spin_system,L,timestep)` |

`evolution` supports all of the following `output` strings. `krylov` supports
all of them except `'total'`, which requires the infinite-time integration
implemented by `evolution`.

| String | Returns |
|---|---|
| `'final'` | final state vector, or a horizontal stack thereof |
| `'trajectory'` | the stack of state vectors over the requested steps |
| `'total'` | integral of the observable trace from start to infinity; requires relaxation and is supported by `evolution` only |
| `'refocus'` | evolves the first vector zero steps, the second one step, the third two, matching indirect evolution after a refocusing pulse |
| `'observable'` | time dynamics of one observable against `coil` |
| `'multichannel'` | several observables as rows, `coil` carrying one column each; a channels-by-time-by-states array for a stack of initial states |

For explicit-propagator Liouville-space GPU `refocus`, the state stack is held in full GPU storage and propagated in batches. Column n still receives n-1 steps. Allow GPU memory for the dense stack and multiplication temporaries, rather than only one state vector.

The optional `destination` argument to `evolution` is a state for
destination-state screening. `krylov` is for the case where `exp(L)` will not
fit in memory but `L` will.

`step` normally computes the matrix-exponential action without forming it. Small
Hilbert density matrices use an explicit exponential; sparse GPU generators are
converted to full GPU storage locally for MATLAB's `expm` compatibility. CPU
generators and the larger-matrix route retain their existing storage behaviour. Pass
one matrix for the centre-point piecewise-constant rule, `{left,right}` for
piecewise-linear, `{left,midpoint,right}` for piecewise-quadratic. Ideal hard
pulses are `step` calls with a pulse operator and the flip angle in place of
the time step:

```matlab
Ly=operator(spin_system,'Ly',parameters.spins{1});
rho=step(spin_system,Ly,rho,pi/2);
```

Coherence selection is analytical rather than by phase cycling:

- `rho=coherence(spin_system,rho,spec)` keeps chosen coherence orders;
  `{{'13C',[1 -1]},{'1H',-1}}` keeps states with ((1 OR -1 on 13C) AND
  (-1 on 1H)). `'electrons'`, `'nuclei'`, and `'all'` work in place of an
  isotope name.
- `rho=correlation(spin_system,rho,orders,spins)` keeps chosen spin
  correlation orders; `orders` is a row vector.
- `rho=homospoil(spin_system,rho,zqc_flag)` emulates a strong homospoil
  gradient, keeping only states of zero frequency with respect to the
  carriers. `zqc_flag` is `'keep'`, which lets zero-quantum coherences
  survive as they do experimentally, or `'destroy'`.
- `[L,rho]=decouple(spin_system,L,rho,spins)` wipes interactions and
  populations involving the named spins until `L` is rebuilt.
- `rho=spinlock(spin_system,Lx,Ly,rho,direction)` applies an analytical
  spin lock.

The first three also run in `zeeman-liouv` and `zeeman-hilb` through an
internal Liouville round trip (`homospoil` there keeps only the
density-matrix diagonal and ignores `zqc_flag`). `coherence` and
`correlation` support Fokker-Planck direct products in every formalism;
`homospoil` supports them in `sphten-liouv` only, because its
`zeeman-liouv` branch takes `sqrt(size(rho,1))` as the density matrix
dimension instead of factoring the spatial part out, which throws for
most spatial dimensions and masks the wrong elements when the combined
row count happens to be a perfect square.

`dnp_freq_scan` honours the polyadic enable switch for `fp-gmres`, using FFT phase differentiation and a spin-sized Fourier-block preconditioner. Nonconverged implicit solves fail explicitly. `fp-backs` retains direct matrix solves, and LvN methods are unchanged. This is a CPU steady-state route, not GPU GMRES.

## The `state` and `operator` grammar

```matlab
rho=state(spin_system,states,spins,method)
A=operator(spin_system,operators,spins,operator_type,format)
```

Both take the same descriptor grammar, in four forms.

| Form | Example | Result |
|---|---|---|
| string, string | `state(spin_system,'Lz','13C')` | sum of the single-spin objects over all spins of that isotope |
| string, vector | `operator(spin_system,'Lz',[1 2 4])` | sum over the numbered spins |
| cell, cell | `state(spin_system,{'Lz','L+'},{1,2})` | the product state LzS+ - density matrix in Hilbert space, state vector in Liouville space |
| projection numbers | `state(spin_system,[-1/2 1/2 0])` | wavefunction formalism only, in a `{'1H','1H','14N'}` system; no third argument |

Valid labels are `'E'` (identity), `'Lz'`, `'Lx'`, `'Ly'`, `'L+'`, `'L-'`,
`'Tl,m'` (irreducible spherical tensor with integer `l` and `m`, e.g.
`'T2,0'`), and the Zeeman-basis central-transition operators `'CTx'`,
`'CTy'`, `'CTz'`, `'CT+'`, `'CT-'`. Valid spin specifications are standard
isotope names plus `'electrons'`, `'nuclei'`, and `'all'`.

`method` in `state` is `'exact'` by default; `'chem'` is a deprecated alias.
Density matrices and Liouville vectors are weighted by `inter.chem.concs` in
both cases; wavefunctions remain unweighted. Detection uses `coil_state` with
all three arguments and no method selector. `operator_type` applies in Liouville
space and is ignored
in Hilbert space: `'comm'` (commutation, the default), `'left'`, `'right'`,
`'acomm'`. `format` is `'csc'` for a MATLAB sparse matrix (default) or
`'xyz'` for a `[rows, cols, vals]` array; adding `'op_cache'` to `sys.enable`
turns on operator caching.

A product of two commutation superoperators is not the commutation
superoperator of the product. In Liouville space you cannot build single-spin
superoperators and multiply them together; ask `operator` for the product
form directly.

Typical usage inside a sequence, and the pattern used throughout
`experiments/`:

```matlab
parameters.rho0=state(spin_system,'Lz',parameters.spins{2});
parameters.coil=coil_state(spin_system,'L+',parameters.spins{2});
Lx_F1=operator(spin_system,'Lx',parameters.spins{1});
Ly_F2=operator(spin_system,'Ly',parameters.spins{2});
```

Thermal equilibrium comes from
`rho=equilibrium(spin_system,I,Q,euler_angles)`, isotropic if `Q` and the
angles are omitted, or from `parameters.needs={'rho_eq'}` and its
context-specific relatives.

## The experiments library

Ready-made sequences live in `experiments/`. The standard signature is
`fid=NAME(spin_system,parameters,H,R,K)`; imaging and spatially encoded
sequences add `,G,F`. A handful take only `(spin_system,parameters)` because
they build their own generators.

**1D and general** - `acquire` (forward evolution and acquisition),
`hp_acquire` (hard pulse), `sp_acquire` (Fokker-Planck soft pulse),
`holeburn`, `cpmg`, `inv_rec`, `sat_rec`, `traject`, `slowpass`
(frequency-domain slow passage), `respiration`, `cp_acquire_soft`,
`cp_contact_hard`, `cp_contact_soft`, `relaxan`, `eqmag`, `fieldsweep`,
`fieldscan_enlev`, `fieldscan_magn`, `rapidscan`.

`slowpass` takes a two-element `sweep` interval in Hz, `npoints>=2`,
`rho0`, and `coil`, and returns the FFT-normalised frequency-domain response.
Relaxation is required and must not be thermalised. Liouville identity
components are excluded only when their sector is decoupled from spin order
in both directions within `tols.liouv_zero`; that sector is shifted to lift
stationary poles, while the spin-order block is unchanged. Spatial contexts
supply `spc_dim` for identity embedding. Selective-reaction coupling retains
the original resolvent; wavefunction inputs bypass identity removal. This
does not regularise other undamped or coupled stationary modes.

`acquire` takes `sweep` (Hz), `npoints`, `rho0`, `coil`, `decouple` (e.g.
`{'15N','13C'}`), and optionally `homodec_oper` with `homodec_pwr` (Hz) and
`dead_time` (seconds). `hp_acquire` adds `pulse_op` and `pulse_angle`
(radians), with optional `echo_time`, `echo_oper`, `echo_angle`.
`sp_acquire` adds `pulse_frq` (Hz, relative to the current rotating frame),
`pulse_phi` (rad), `pulse_pwr` (rad/s), `pulse_dur` (s), `pulse_rnk`
(Fokker-Planck cut-off rank - start at 2 and increase until the answer stops
moving), and `method` (`'expv'`, `'expm'`, `'evolution'`). When `sim2liouv`
admits a Hilbert-space sequence, it converts `homodec_oper` into a
commutation superoperator alongside the other operator-like fields.

**Liquid-state 2D** (`experiments/nmr_liquids/`) - `cosy`, `gcosy`,
`ct_cosy`, `dqf_cosy`, `ecosy`, `tocsy`, `noesy`, `roesy`, `hoesy`,
`noesyhsqc`, `hsqc`, `clip_hsqc`, `ct_hsqc`, `hmqc`, `hmbc`, `hetcor`,
`coloc`, `inept`, `dept`, `deptq`, `inadequate`, `inadequate_2d`, `mqs`,
`mqs_refocus`, `crazed`, `pansy_cosy`, `pansy_triple`.

`hsqc` takes `sweep` and `npoints` as `[F1 F2]` vectors, `spins` as `{F1 F2}`
nuclei, `parameters.J` (working scalar coupling, Hz), `decouple_f2` (nuclei
decoupled in F2) and `decouple_f1` (nuclei receiving midpoint 180-degree
refocusing pulses in F1); it defaults `rho0` to `'Lz'` and `coil` to `'L+'`
on `parameters.spins{2}`. `noesy` takes `parameters.tmix` (mixing time,
seconds) and `parameters.oldschool` to disable the homospoil gradient before
mixing. `gcosy` takes `parameters.angle` (second pulse angle in radians,
allowing COSY45 and COSY60), `g_amp` (Gauss/cm, default 3), `g_dur` (seconds,
default 2e-3), `g_stab_del` (default 2e-4), `s_len` (active sample length in
cm, default 1.5), and `parameters.pathway`, one of `'P'` (default), `'N'`, or
`'P+N'`.

Natural-abundance simulations should use isotope dilution -
`subsystems=dilute(spin_system,isotope,tuples)` returns a cell array of
spin systems, one per isotopomer, to be looped over. `kill_spin` (also used by `dilute`) drops compiled records touching a removed particle and reindexes surviving coupling tensors, proximity, relaxation, and bosonic channels, including nested spin-modulation targets. It clears mode assumption strengths and removes the mode container when no C, V, or T particles remain. An existing basis is rebuilt; if there was no basis, call `basis` explicitly. Reapply `assume` before constructing a Hamiltonian. After removing the final mode, ordinary spin-only assumption sets and retention options apply.

**Bruker ports** (`experiments/bruker/`) - literal translations of the
echo/antiecho gradient-selected pulse programs: `hmqcetgp`, `hmqcetgpsi`,
`hsqcetgp`, `hsqcetgpsi`, `hsqcedetgp` (multiplicity-edited); the `si`
variants are sensitivity-improved.

**Protein** (`experiments/nmr_protein/`) - `hnca`, `hnco`, `hncoca`,
`hncaco`, `hcanh`, `hcch_cosy`, all with pre-set protein J-couplings.

**Solids** (`experiments/nmr_solids/`) - `cn2d_dq`, `cn2d_sq`, `dante`,
`fslghetcor`, `mqmas`, `pdsd`, `redor`, `wise`. Overtone and NQR variants are
in `experiments/overtone/` (`overtone_a`, `overtone_cp`, `overtone_dante`,
`overtone_pa`) and `experiments/nqr/nqr_pa.m`.

**EPR dipolar** (`experiments/esr_dipolar/`) - `deer_3p_hard_deer`,
`deer_3p_hard_echo`, `deer_3p_soft_deer`, `deer_3p_soft_hole`,
`deer_4p_soft_deer`, `deer_4p_soft_hole`, the diagnostic suites
`deer_3p_soft_diag(spin_system,parameters)` and `deer_4p_soft_diag`, the
analytical two-spin trace `deer=deer_analyt(D,J,t)`, plus `eseem`,
`oopeseem`, `ridme`, and `sifter`.

**EPR hyperfine** (`experiments/esr_hyperfine/`) - `endor_cw` (fast
isotropic CW-ENDOR), `endor_davies` (explicit soft pulses, orientation
selection), `endor_mims`, `endor_mims_ideal`, `endor_mims_echo`, `hyscore`.
In `endor_mims_ideal`, the spectrum has the same shape as the RF offset grid
`parameters.n_frq`.

**Hyperpolarisation** (`experiments/hyperpol/`) - `dnp_field_scan`,
`dnp_freq_scan`, `dnp_time_dep`, `beamdnp`, and the pairs `noveldnp`,
`topdnp`, `xixdnp` with their `_steady` steady-state versions.
`solid_effect(spin_system,parameters)` and `masdnp(spin_system,parameters)`
build their own generators and take no `H,R,K`.

**Imaging** (`experiments/imaging/`, signature `(…,H,R,K,G,F)`) -
`basic_1d_hard`, `slice_select_1d`, `spin_echo`, `grad_echo`,
`phase_enc_2d`, `phase_enc_3d`, `epi_2d`, `epi_3d`, `fse`, `spiral`,
`cpmg_dec`, `udd_dec`, `dpfgse_select`, `dpfgse_suppress`, `press_1d`,
`press_2d`, `press_voxel_1d`, `press_voxel_2d`, `press_voxel_3d`. The
parameter idiom, from `phase_enc_2d`: `t_echo` (seconds), `pe_grad_amp` and
`ro_grad_amp` (T/m), `pe_grad_dur` and `ro_grad_dur` (seconds), `image_size`,
and the optional `diff_g_amp`.

**Spatially encoded and diffusion** (`experiments/spen/`, also
`(…,H,R,K,G,F)`) - `psyche`, `psycosy`, `spencosy`, `spendosy`,
`spendosycosy`, `ufmq`, `dosy_oneshot`, `idosyzs`, `st_ideal`.

**Other** - `experiments/zulf/` (`zerofield`, `zulf_abrupt`),
`experiments/singlets/` (`m2s`, `s2m`, both
`(spin_system,L,Hx,Hy,rho,J,delta_v)`), `experiments/spin_chem/` (`rydmr`,
`rydmr_exp`), `experiments/microfluidics/simple_flow` (meshflow context),
`experiments/rdc/`, and `experiments/pseudocon/`.

## Acquisition and processing conventions

For nD experiments `parameters.sweep`, `parameters.npoints`, and
`parameters.zerofill` become vectors ordered `[F1 F2]`, and
`parameters.spins` lists the channels in the same order. The FID array is
laid out with dimension 1 as F2 (direct) and dimension 2 as F1 (indirect),
so `zerofill(2)` is used for the F2 transform and `zerofill(1)` for F1.

Quadrature in the indirect dimension is arithmetic in the calling script,
not a library function. Sequences that produce States data return `fid.cos`
and `fid.sin`; sequences that produce echo/anti-echo data return `fid.pos`
and `fid.neg`.

```matlab
% States (hypercomplex)
fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}});
fid.sin=apodisation(spin_system,fid.sin,{{'sqcos'},{'sqcos'}});
f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1));
f1_sin=real(fftshift(fft(fid.sin,parameters.zerofill(2),1),1));
f1_states=f1_cos-1i*f1_sin;
spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2);

% Echo/anti-echo
f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1);
f1_neg=fftshift(fft(fid.neg,parameters.zerofill(2),1),1);
fid=f1_pos+conj(f1_neg);
spectrum=fftshift(fft(fid,parameters.zerofill(1),2),2);
```

Magnitude-mode 2D data are transformed in one call, with the argument order
`fft2(fid,parameters.zerofill(2),parameters.zerofill(1))`.

### Apodisation

```matlab
fid=apodisation(spin_system,fid,winfuns,fp_half)
```

`winfuns` is a cell array of cell arrays, one entry per FID dimension,
indexed by absolute dimension number; an empty cell `{}` means "do nothing in
this dimension" and also excludes it from first-point halving. `fp_half`
defaults to `true()`; set it to `false()` when applying several windows to
the same dimension in sequence, so the first point is halved only once.

| Window | Definition |
|---|---|
| `'none'` | ones; the first point is still halved |
| `'crisp'` | cos(x)^8 half-bell over [0,pi/2] |
| `'cos'` | cos(x) over [0,pi/2] |
| `'sqcos'` | cos(x)^2 over [0,pi/2] |
| `'sin'` | sin(x) over [0,pi] |
| `'sqsin'` | sin(x)^2 over [0,pi] |
| `'exp',k` | exp(-k*x), x in [0,1]; k e-folds of decay by the last FID point |
| `'gauss',k` | exp(-k*x^2), x in [0,1] |
| `'kaiser',k` | `kaiser(npts,k)`, sidelobe attenuation |
| `'bad-z1',k` | misset Z1 shim (sinc); k=10 is a reasonable 1H 600 MHz guess |
| `'bad-z2',k` | misset Z2 shim (Fresnel); k=40 is a reasonable guess |

The five parameterised windows must be given as `{'name',k}` with `k` a
finite real scalar; the other six error if handed a parameter.

### Fourier transform

There is no FT wrapper - MATLAB's `fft`, `fft2`, and `fftshift` are called
directly. The 1D idiom is

```matlab
fid=apodisation(spin_system,fid,{{'exp',6}});
spectrum=fftshift(fft(fid,parameters.zerofill));
```

Explicit axes, when needed: `ax=ft_axis(offset,sweep,npoints)`,
`[f_shift,f,df]=fft_freq_axis(npts,dt,zf)`,
`[t_shift,t,dt]=ifft_time_axis(npts,df,zf)`,
`axis_hz=sweep2ticks(offs,sweep,npoints)`.

### Plotting

```matlab
plot_1d(spin_system,spectrum,parameters,varargin)
[axis_f1,axis_f2,spectrum]=plot_2d(spin_system,spectrum,parameters,...
                                   ncont,delta,k,ncol,m,signs)
```

`plot_1d` passes `varargin` through to MATLAB's `plot`, splits complex
spectra into real and imaginary traces with a legend, and requires
`size(spectrum,1)` to equal `parameters.zerofill` (filled from
`parameters.npoints` if `zerofill` is absent). It reads `parameters.sweep`,
`parameters.offset`, `parameters.spins` (exactly one element),
`parameters.axis_units`, `parameters.derivative`, and
`parameters.invert_axis`.

Axis units are converted in `axis_1d`, which accepts `'ppm'`, `'Gauss'`,
`'mT'`, `'Hz'`, `'kHz'`, `'MHz'`, `'MHz-labframe'`, `'GHz'`,
`'GHz-labframe'`, `'gtensor'`, and `'points'`. Magnetic-field units use the
free-electron g-tensor, and `'ppm'` errors if the magnet field is zero.
`plot_2d` accepts a shorter list - `'ppm'`, `'Gauss'`, `'Hz'`, `'kHz'`,
`'MHz'`, `'points'` - and `plot_3d` only `'ppm'`, `'Gauss'`, `'Hz'`.

All nine `plot_2d` arguments are mandatory. `delta` is a four-vector
`[posmin posmax negmin negmax]` in [0,1], ascending within each pair;
`ncont`, `k`, `ncol`, and `m` are positive integers; `signs` is
`'positive'`, `'negative'`, or `'both'`. A working call:

```matlab
plot_2d(spin_system,real(spectrum),parameters,20,[0.02 0.2 0.02 0.2],2,256,6,'both');
```

`plot_2d` transposes internally and reverses both axes, so the F2-in-dimension-1
convention above comes out the right way round. In `slice_2d`, clicking a point selects a row for the F1 trace and a column for the F2 trace. Related plotters:
`plot_3d(spin_system,spectrum,parameters,nsurf,delta,k,signs)`,
`stack_2d(spin_system,spectrum,parameters,stack_dim,alpha_fun)`, `slice_2d`
and `int_2d` (as `plot_2d`, `int_2d` with a trailing filename), `plot_uf`,
and `mri_2d_plot(mri,parameters,method)` with `method` in
`{'image','phantom','k-space'}`. Cropping is
`[spec,parameters]=crop_2d(spin_system,spec,parameters,crop_ranges)` with
`crop_ranges={[f1_min f1_max],[f2_min f2_max]}` in ppm only. Figures use the
house wrappers `kfigure`, `ktitle`, `kgrid`, `kxlabel`, `kylabel`, `kzlabel`,
`klegend`, `kletter`, `kcolourbar`, and `scale_figure`. `kletter(letter_label)`
takes a one-element character string and draws a bold panel letter inside the plot
box, with its left edge and its cap line 10 points from the left and the top
edge of the box; the offsets are computed when it is called, so call it after the
figure has its final size and a tiled layout all of its tiles; `fig2tiles` re-applies the tagged letters on the
retiled axes of a merged figure.

## Cluster-expansion context

`[answer,diagnostics]=clexpn(spin_system,pulse_sequence,parameters,assumptions)` accepts `create` output without a full-bath basis. Supply ordered physical `parameters.observed` IDs and strictly increasing, nonnegative total `parameters.times` in seconds. The registered experiment receives a retained `(cluster_system,parameters,H,R,K)`; all observed spins are retained in every cluster. The first output is the numerical experiment answer, not an engine record.

Registered handles are `@free_evolution`, `@hahn_echo`, `@cpmg_echo`, `@bath_pumped_hahn`, `@cce_eseem`, `@finite_drive`, `@refocusing_train`, and `@noise_correlation` in `experiments/cce/`. `cce_capabilities` declares admissible generators, output and reduction algebra; arbitrary callbacks are not accepted just because their signatures match. The context builds retained systems and generators; experiments own pulse timing and propagation.

Physical tensors remain canonical `create` data. Preparation, detector, transition, pulse, drive, pump, and jump fields are experiment inputs. Cluster approximation, bath IDs, graph, order, partition, sampled states, numerical policy, scheduling, and checkpoints belong under `parameters.cce`. `cce_options` rejects unknown fields and incompatible ownership. See [Cluster correlation expansions](cce.md) for the complete native guide and supported output restrictions.

The default approximation is `generalized`. Projected calculations require an explicit observed transition; arbitrary finite drives and observed-density output require generalised evolution. Dissipative capability is registered for free evolution, Hahn echo, and finite drive, not every experiment. Native `bath_jumps` declare ordered support, operator, and rate; conventional dissipation uses retained Liouville space, while projected branches use an explicit conditional representation. Exchange is not supported. The observation compiler returns the actual initial detector reference, and the experiment applies normalisation once.

For a projected coherent `@free_evolution` signal, an explicit `parameters.coil` must be proportional to the selected transition detector. Its reported floating-point propagation estimate includes the receiver magnitude and, for normalised signals, division by the magnitude of the initial signal. A valid full-bath result need not meet the requested numerical budget.

`diagnostics` contains validity, numerical and sampling uncertainty, conditioning, policy, convergence, and execution evidence. Invalid values are not repaired with arbitrary denominator floors. The `cce.averaging` setting selects reduction before ensemble averaging or raw-record averaging before reduction; additive noise correlation and multiplicative coherence use distinct declared algebras. Omitted cluster order includes the selected bath, but does not establish convergence in bath extent, graph selection, or state sampling.

### Orientation grids

For one orientation, set `parameters.orientation` to a radian Euler triple. For a quadrature, supply `parameters.grid` as a shipped grid name or an explicit structure containing `alphas`, `betas`, `gammas` and `weights`. Call `clexpn(spin_system,pulse_sequence,parameters,assumptions)` directly in both cases. Set `parameters.sum_up=false` to return a column cell array of unweighted responses in input order. Zero-weight nodes are skipped; `diagnostics.grid_indices` maps each response to its original grid entry. `parameters.cce.averaging` distinguishes reducing each orientation before averaging from averaging matching raw cluster records before reduction.

Each orientation runs independent complete dynamics. The task's cce.averaging selects expansion-first averaging or raw-record averaging before expansion with a common cluster family; the observed-only reference participates in both. Hamiltonians are never averaged. Orientations are serial; nested cluster/configuration scheduling may use the existing pool, but becomes serial inside an outer worker.

The clexpn native boundary rejects parameters.sequence, including empty arrays. Use finite_drive, refocusing_train, or a named echo experiment with native physical descriptors; free_evolution is uncontrolled.
