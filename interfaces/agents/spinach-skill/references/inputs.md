# Spinach input specification

## Contents

- [The `sys` structure](#the-sys-structure)
- [Isotope naming](#isotope-naming)
- [Zeeman interactions](#zeeman-interactions)
- [Couplings](#couplings)
- [Coordinates, dipolar couplings, periodic boundaries](#coordinates-dipolar-couplings-periodic-boundaries)
- [Magnetic susceptibility and paramagnetic shifts](#magnetic-susceptibility-and-paramagnetic-shifts)
- [Giant spin (ligand field) terms](#giant-spin-ligand-field-terms)
- [Chemistry and kinetics](#chemistry-and-kinetics)
- [Partial ordering, temperature and relaxation inputs](#partial-ordering-temperature-and-relaxation-inputs)
- [The `bas` structure](#the-bas-structure)
- [Unit conversions on the way in](#unit-conversions-on-the-way-in)
- [Importing from quantum chemistry](#importing-from-quantum-chemistry)
- [Importing structures and databases](#importing-structures-and-databases)

Everything a simulation knows about physics enters through two structures,
`sys` and `inter`, which are consumed by `create`, and one structure, `bas`,
consumed by `basis`. There are no defaults for physical quantities: a missing
interaction is either an explicit warning ("zeros assumed") or an error, never
a guess.

```matlab
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
```

`create` refuses to run when called directly from MATLAB's base workspace; a
named script or function adds a caller stack frame and is accepted. Keep
simulations in a named script or function rather than calling `create`
interactively. Unrecognised fields in `sys` are fatal: `create` strips the
fields it understands and errors on whatever is left.

## The `sys` structure

| Field | Type | Meaning |
|---|---|---|
| `sys.isotopes` | cell array of strings, one per spin | Mandatory. Particle specification, order defines spin numbering. |
| `sys.magnet` | real scalar | Mandatory. Magnetic induction in tesla. Zero is legal (zero-field work). |
| `sys.labels` | cell array of strings, same length as `sys.isotopes` | Optional text labels; printed in diagnostics and resolved to indices by `idxof`. |
| `sys.output` | `'hush'`, `'console'`, or a file name | Destination of the console report. Default is the console. |
| `sys.scratch` | directory path | Scratch directory; defaults to `<root>/scratch`, created if absent. |
| `sys.disable` | cell array of strings | Switches off internal algorithms. |
| `sys.enable` | cell array of strings | Switches on optional algorithms. |
| `sys.tols` | structure | Overrides numerical cut-offs. |
| `sys.parallel` | `{pool_type,nworkers}` | Parallel pool specification, e.g. `{'processes',8}`; defaults to a local process pool leaving one core to the OS. |
| `sys.parprops` | cell array of name-value pairs | Extra properties passed to the parallel cluster object. |

Legal `sys.disable` entries, anything else being an error: `'pt'` (non-interacting subspace detection),
`'symmetry'`, `'krylov'`, `'clean-up'`, `'hygiene'` (start-up health checks),
`'dss'`, `'expv'`, `'trajlevel'`, `'merge'`, `'colorbar'`, `'asyredf'`. Legal
`sys.enable` entries: `'zte'` (zero track elimination, off by default), `'gpu'`, `'op_cache'`, `'ham_cache'`, `'prop_cache'`,
`'greedy'`, `'paranoia'` (tight tolerances), `'cowboy'` (loose
tolerances), `'polyadic'`, `'sodd'` (spin-orbit corrections to dipolar
couplings), `'dafuq'`. With `'polyadic'` enabled, `v2fplanck` supports scalar and voxel-wise velocities; spin-space Kronecker extension preserves prefactors by nesting affixed polyadics.

`sim2liouv` refreshes existing basis cache identities after Hilbert-to-Liouville
conversion; cached operators and Hamiltonians remain representation-specific.

`sys.tols` subfields are listed and defaulted in `tolerances.m`. The two that
change physics rather than performance are `inter_cutoff`, below which coupling
tensors are discarded (2-norm, in Hz), and `prox_cutoff`, which bounds the
distance over which spins count as proximate.

## Isotope naming

`spin(name)` is the database behind `sys.isotopes`; it returns the magnetogyric
ratio in rad/(s·tesla) and the multiplicity.

| Specification | Particle |
|---|---|
| `'1H'`, `'13C'`, `'15N'`, `'195Pt'`, ... | Nucleus: mass number followed by element symbol. |
| `'E'`, `'E+'` | Electron and positron, multiplicity 2. |
| `'E4'`, `'E16'` | High-spin electron; the integer is the **multiplicity**, so `'E4'` is S=3/2 and `'E16'` is S=15/2. |
| `'G'` | Ghost spin, gamma=0 and multiplicity 1; a placeholder that carries coordinates but no magnetism. |
| `'N'`, `'antiN'`, `'M'`, `'M+'` | Neutron, antineutron, negative and positive muons. |
| `'anti1H'` | Antiproton; `'1H'` remains the proton row. |
| `'Lambda'`, `'Sigma+'`, `'Sigma-'`, `'Xi0'`, `'Xi-'`, `'Omega-'` | Static-moment hyperons; `anti`-prefixed keys name their CPT-qualified counterparts. |
| `'99Tc_m'`, `'180Ta_m'`, ... | Distinct evaluated nuclear states, not unqualified ground-state aliases. |
| `'C#'`, `'V#'`, `'T#'` | Cavity mode, phonon mode, transmon; the integer is the number of levels, and gamma is zero for all three. |

The fixed signature is `[gamma,multiplicity,data]=spin(name)`. Existing one-/two-output
calls retain their numerical API. A third output returns a sourced metadata row;
`[~,~,isotopes]=spin('table')` exposes all physical rows, including unknown or
tentative properties that are not usable as confirmed simulation inputs. Missing
moments/spins raise a data-unavailable error, never fabricated zeros. Abundance is
a fraction, quadrupole moment is in barns, and half-life is in seconds; metadata
does not automatically apply abundance weights or radioactive decay. Interval-only
abundance has no invented midpoint. The editable literature TSV and offline MAT
builder are documented in `etc/isotopes_sources.md`; the runtime table loads once
per process and `clear spin` reloads an updated payload. Natural tantalum-180 is the
`180Ta_m` isomer, not the short-lived `180Ta` ground state. Antiparticle values under
CPT are explicitly qualified, not independent measurements. Spin below one forbids
a spectroscopic quadrupole moment; higher-spin missing values, including Omega,
remain unknown rather than zero. Nuclear states with adopted half-lives below one
second (values, estimates, or upper limits) are excluded; stable nuclei, exactly
one second, unknown lifetimes, unresolved lower limits, and magnetic particles
are retained. Alternative environmental lifetimes do not override this cutoff.

`iselectron` and `isnucleus` test a specification string. `isoswap(sys,inter,
spins,new_iso)` performs isotope replacement and rescales all interactions
accordingly, wiping quadratic and higher-order couplings with a warning.

`isot2elem({'1H','13C','35Cl'})` returns `{'H','C','Cl'}`, retaining the input
cell-array shape and order. It strips digits only, without an isotope lookup;
isomer suffixes and other non-digit characters are not removed.

## Zeeman interactions

Three mutually compatible specifications exist; whatever is supplied is summed
into a single tensor per spin. Nuclear values are chemical shifts in **ppm**;
electron values are **g-tensors in Bohr magneton units**. Getting this
distinction wrong is silent, because both are dimensionless numbers.

```matlab
inter.zeeman.scalar={1.0 1.5};                 % isotropic, 1 x nspins cell
inter.zeeman.eigs={[10 20 30] []};             % principal values, 1 x 3 each
inter.zeeman.euler={[0 pi/4 0] []};            % ZYZ active, radians
inter.zeeman.matrix={eye(3)*2.0 []};           % full 3x3 tensors
```

`scalar` has one element per spin; empty or zero entries are skipped. `eigs`
and `euler` must be cell arrays of identical size with **identical non-empty
patterns** — supplying one without the other, or leaving `euler` empty where
`eigs` is not, is an error; use `[0 0 0]` for no rotation. `matrix` must be a
`1 x nspins` cell array of real `3x3` matrices or empties; a column cell array
is rejected.

Internally, nuclear tensors become `(eye(3)+1e-6*matrix)*basefrq` and electron
tensors become `matrix*basefrq/g_free`, where `basefrq=-gamma*sys.magnet`. If
neither `inter.zeeman` nor `inter.suscept` is present, bare magnet frequencies
are assumed and a warning is printed.

## Couplings

All spin-spin couplings — J, dipolar, hyperfine, quadrupolar, zero-field
splitting — live in `inter.coupling` and are all in **hertz** on input.
The three specifications superpose, so a J-coupling and a dipolar tensor for
the same pair simply add.

```matlab
inter.coupling.scalar=cell(nspins,nspins);     % Hz, isotropic
inter.coupling.scalar{1,2}=7.0;

inter.coupling.eigs=cell(nspins,nspins);       % Hz, principal values
inter.coupling.euler=cell(nspins,nspins);      % radians, ZYZ active
inter.coupling.eigs{1,2}=[-1e3 -1e3 2e3];
inter.coupling.euler{1,2}=[0 0 0];

inter.coupling.matrix=cell(nspins,nspins);     % Hz, full 3x3 tensors
inter.coupling.matrix{1,2}=A;
```

All three are `nspins x nspins` cell arrays. `eigs` and `euler` must have the
same dimensions and identical non-empty patterns. Only one triangle needs
filling for a given pair; the kernel handles the rest. Couplings whose norm is
below `sys.tols.inter_cutoff` are dropped, as are all couplings involving ghost
spins. Couplings between spins that belong to different chemical subsystems are
a hard error.

**Diagonal elements are single-spin quadratic couplings.** `{n,n}` holds the
quadrupolar tensor of nucleus *n* or the zero-field splitting tensor of
electron *n*:

```matlab
inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0]);   % C_q, eta, I, eulers
inter.coupling.matrix{2,2}=zfs2mat(D,E,alp,bet,gam);          % D, E in Hz
```

`inter.ignore` is a cell array of two-element index vectors; each listed pair
has its coupling tensor deleted after all other processing.

## Coordinates, dipolar couplings, periodic boundaries

```matlab
inter.coordinates={[0.00 0.00 0.00]
                   [0.00 0.00 1.00]};
```

A cell array of real `1x3` vectors in **angstrom**, one per spin, empties
allowed. Supplying it invokes the dipolar module, which computes every point
dipolar coupling from geometry — so coordinates and explicit dipolar tensors
for the same pair will double-count. Without coordinates, a warning is printed
and point dipolar interactions are taken as zero.

`inter.pbc` is a cell array of lattice translation row vectors in angstrom;
supplying it makes the dipolar module apply periodic boundary conditions with
lattice summation. An empty cell array means a standalone system.
`cubic_lattice(isotope,spacing,n_periods)` returns `sys` and `inter` with
`isotopes`, `coordinates` and `pbc` already set.

To compute a tensor by hand instead: `[d,alp,bet,gam,M]=xyz2dd(r1,r2,isotope1,
isotope2)` returns the dipolar coupling constant and tensor in **rad/s** from
coordinates in angstrom, using free-particle magnetogyric ratios;
`A=xyz2hfc(exyz,nxyz,isotope)` returns a point-dipole hyperfine tensor in
**gauss** and is the one to use when an electron is involved. `conmat(xyz,r0)`
builds a connectivity matrix; `nearest_spin` and `dihedral` are geometry
helpers.

## Magnetic susceptibility and paramagnetic shifts

```matlab
inter.suscept.chi={[ 0.0883 -0.0904  0.0822
                    -0.0904 -0.1011 -0.0149
                     0.0822 -0.0149  0.0128]};
inter.suscept.xyz={[10.0 2.5 3.9]};
```

`chi` is a cell array of `3x3` susceptibility tensors in **cubic angstrom**;
`xyz` gives the position of each centre in angstrom. For every nucleus that has
coordinates, `create` computes the paramagnetic shielding tensor with `xyz2pms`
and adds it, in ppm, to that nucleus's Zeeman tensor; electrons are unaffected.
Quantum chemistry quotes susceptibility in cgs-ppm (cm³/mol), so convert with
`cgsppm2ang` first. `chi=g2chi(g,T,S)` is the high-temperature Curie estimate
from a g-tensor, `T` in kelvin.

## Giant spin (ligand field) terms

Available for electron spins only. `inter.giant.coeff{n}{k}` is the vector of
`2k+1` spherical-tensor coefficients of rank `k` for spin `n`, in hertz;
`inter.giant.euler{n}{k}` is the corresponding `1x3` Euler angle vector in
radians. Both cell arrays must have one entry per spin, and every rank present
in `coeff` must have its Euler angles supplied. With `ham_cache` enabled,
giant-spin coefficients and retention strengths distinguish cache entries,
including full and Zeeman-only Hamiltonians.

```matlab
inter.giant.coeff={{[0 0 0],Bkq{2},[0 0 0 0 0 0 0],Bkq{4}}};
inter.giant.euler={{[0 0 0],[0 0 0],[0 0 0],[0 0 0]}};
```

Stevens operator coefficients are converted to the irreducible spherical tensor
convention with `stev2sph(k,Bkq)`.

## Chemistry and kinetics

```matlab
inter.chem.parts={1,2};                        % one spin per species
inter.chem.concs=[1.0 1.0];
inter.chem.reactions={struct('reactants',1,'products',2,...
                            'matching',[1 2],'rate',2e4),...
                      struct('reactants',2,'products',1,...
                            'matching',[2 1],'rate',2e4)};
```

- `parts` — cell array of numeric index vectors, one per chemical subsystem,
  disjoint and within the spin count; an empty entry denotes a spin-free substance.
  Reaction-bearing inputs require explicit parts; otherwise the default is one
  subsystem containing everything.
- `concs` — non-negative initial concentrations, one per subsystem; required
  for multiple substances. Propagated unit coordinates carry subsequent
  concentrations. `state` weights each block by its initial concentration;
  `coil_state(...,'exact')` constructs unweighted detection vectors.
- `reactions` — cell array of scalar records with row-vector `reactants` and
  `products`, two-column global-spin `matching`, and non-negative scalar or
  time-dependent `rate`. First-order rates are in inverse seconds; order-m rates
  use concentration^(1-m)/s. Matched spins must have identical isotopes. Empty
  products denote untracked loss; repeated substance indices give stoichiometry.
  Unmatched source spins are traced out and new product spins arrive unpolarised.
- `closure` — per-record `additive` default, or `product` to retain cross-reactant
  polarisation products. Higher-order rates multiply the other reactant
  concentrations in the instantaneous state; their units depend on reaction order.
- `selector` — optional named singlet/triplet channel (including Jones–Hore
  variants) and two electron indices on a single reactant, or a pair of local
  left/right projector matrices. See the selective-loss recipe in relaxation.md.

Spin replacement uses explicit matching records, not a flux matrix.
Legacy rates, flux, and radical-pair fields are rejected. Reaction records are supported in both Liouville formalisms. `zeeman-hilb` supports first-order matrix actions; `zeeman-wavef` is storage-only.
`merge_inp(sys_parts,inter_parts)` combines `sys`/`inter` structures from
separate DFT calculations into one input set, offsetting spin and subsystem
indices; non-extensive fields such as `magnet` and `temperature` must agree
across parts.

## Partial ordering, temperature and relaxation inputs

`inter.order_matrix` is a cell array of `3x3` Saupe order matrices, one per
chemical subsystem, used by `residual` for RDC and RACS work in liquid
crystals. `inter.temperature` is in kelvin; if absent, 298 K is assumed with a
warning.

Relaxation enters through `inter.relaxation`, `inter.rlx_keep`,
`inter.equilibrium`, `inter.tau_c`, `inter.r1_rates`, `inter.r2_rates`,
`inter.damp_rate`, `inter.lind_*`, `inter.weiz_*`, `inter.nott_*`,
`inter.srfk_*`, `inter.srsk_sources` and `inter.rlx_dfs`, all covered in
`references/relaxation.md`. Rates are in hertz, correlation times in seconds.

## The `bas` structure

```matlab
bas.formalism='sphten-liouv';
bas.approximation={'IK-2'};
bas.connectivity={'scalar_couplings'};
bas.prox_level={3};
```

Every field except `formalism` is a cell with exactly one entry per chemical
substance, even when there is only one. The table gives the contents of each
entry (the filter rows already describe the outer cell). There is no scalar
broadcast. `manual{n}` has local spin columns; `sym_spins{n}` is a cell of local
spin-index vectors. Numeric longitudinal and zero-quantum filter labels remain
global. Empty depth/connectivity entries are used where the local approximation
does not need that setting. The compiled descriptors are `bas.basis{n}`, with
unit rows first and `bas.offsets` delimiting the direct-sum blocks.
Coherence, correlation, homospoil, and decoupling selections retain global
spin labels at their interface; they map those labels into each substance’s
local descriptor and apply the resulting masks at its offset.

| Field | Legal values | Notes |
|---|---|---|
| `formalism` | `'sphten-liouv'`, `'zeeman-liouv'`, `'zeeman-hilb'`, `'zeeman-wavef'` | Mandatory. |
| `approximation` | `'none'`, `'IK-0'`, `'IK-1'`, `'IK-2'`, `'IK-DNP'`, `'IK-SBS'` | Mandatory. Only `'none'` is legal outside `sphten-liouv`. |
| `connectivity` | `'scalar_couplings'`, `'full_tensors'` | Required by, and only legal for, `IK-1`, `IK-2`, and `IK-SBS`. `IK-1` and `IK-2` are spin-only and refuse systems with bosonic modes. In `IK-SBS`, bosonic mode couplings above `tols.inter_cutoff` (pairwise channels, and the spin pairs and spins modulated through `inter.modes.coupling_mod` and `inter.modes.zeeman_mod`, linked to their modes) are added to the coupling graph. |
| `inter_level` | positive integer; `1x3` integer vector for `IK-DNP` and `IK-SBS` | Required by `IK-0`, `IK-1`, `IK-DNP`, `IK-SBS`. Cannot exceed the number of spins; clipped to the spin count of each chemical substance. For `IK-DNP` the three entries bound electrons, spins and nuclei respectively. For `IK-SBS` they are the correlation levels on the boson-boson, spin-boson, and spin-spin coupling graphs; inside spin-boson subgraphs, pure boson-boson correlations above the first level and pure spin-spin correlations above the third level are dropped. |
| `prox_level` | positive integer | Required by, and only legal for, `IK-1` and `IK-2`. Clipped to the spin count of each chemical substance. |
| `projections` | cell array with one row vector of integers per chemical substance | Keeps only the listed total projection quantum numbers in that substance; an empty element means no filter. `sphten-liouv` only. Single substance: `bas.projections={+1}`. |
| `longitudinal` | cell array with one cell array of isotope strings or spin index vectors per chemical substance | Keeps only longitudinal states on those spins of that substance. `sphten-liouv` only. Single substance: `bas.longitudinal={{'15N'}}`. |
| `zero_quantum` | cell array with one cell array of isotope strings or spin index vectors per chemical substance | Keeps only states that are zero-quantum over the union of the listed spins of that substance. `sphten-liouv` only. Single substance: `bas.zero_quantum={{'1H'}}`. |
| `manual` | logical matrix with `numel(chem.parts{n})` columns | Explicit local subgraph list, one subgraph per row. |
| `sym_group` | cell array from `S2`, `S3`, `S4`, `S4A`, `S5`, `S6`, `S6A`, `S8A` | Permutation symmetry groups. |
| `sym_spins` | cell array of index vectors | One vector per group, at least two spins each, no spin in two groups, no group spanning two chemical substances. Mandatory alongside `sym_group`. |
| `sym_a1g_only` | logical | Keep only the fully symmetric irreducible representation. |

`sphten-liouv` refuses multiplicities above 16. `IK-DNP` requires both
electrons and nuclei and nothing else in the system. `IK-SBS` requires both
spins and bosonic modes (`C`, `V`, or `T` particles).

## Unit conversions on the way in

Spinach works in rad/s internally and converts on absorption, so every helper
below produces a number in the units `create` expects.

| Helper | Converts |
|---|---|
| `gauss2mhz(hfc_gauss,g)` / `mhz2gauss(hfc_mhz,g)` | Hyperfine couplings gauss ↔ MHz. `g` optional, free-electron value by default. |
| `mt2hz(hfc_mt,g)` | Hyperfine couplings millitesla → Hz. |
| `g2freq(g,B)` | g-value and field in tesla → electron Zeeman frequency in Hz. |
| `ppm2hz(ppm,B0,nucleus)` / `hz2ppm(hz,B0,nucleus)` | Chemical shift ↔ offset, signs of magnetogyric ratios preserved. |
| `hz2icm` / `icm2hz` | Hz ↔ cm⁻¹, for parameters quoted spectroscopically (ORCA prints ZFS in cm⁻¹). |
| `eeqq2nqi(C_q,eta_q,I,eulers)` | C_q in Hz and asymmetry → `3x3` quadrupolar tensor in Hz. |
| `castep2nqi(V,Q,I)` | CASTEP EFG in atomic units and quadrupole moment in barn → `3x3` NQI tensor in Hz. |
| `weblab2nqi(C_q,eta_q,I,alpha,theta,phi)` | Weblab one-cone model → NQI tensors in Hz. |
| `zfs2mat(D,E,alp,bet,gam)` | D and E in Hz → symmetric `3x3` matrix in Hz. |
| `anas2mat(iso,an,as,alp,bet,gam)`, `axrh2mat(iso,ax,rh,alp,bet,gam)`, `spsk2mat(iso,sp,sk,alp,bet,gam)` | Haeberlen-Mehring anisotropy/asymmetry, axiality/rhombicity, and Herzfeld-Berger span/skew → `3x3`. `mat2axrh(M)` is the inverse of the second. |
| `cgsppm2ang` / `ang2cgsppm` | Susceptibility cgs-ppm ↔ Å³. |
| `euler2dcm` / `dcm2euler` | ZYZ active Euler angles ↔ direction cosine matrix, radians. |
| `frac2cart(a,b,c,alp,bet,gam,ABC)` | Fractional crystallographic → Cartesian coordinates; cell angles in degrees. |
| `fwhm2rlx(fwhm)` | Line width in Hz → approximate R2 in Hz. |

Two mistakes recur. Hyperfine couplings from EPR literature are usually quoted
in gauss or millitesla and must be converted before they enter
`inter.coupling`; and a coupling supplied as a full `3x3` matrix is still in
hertz, not rad/s, however large the numbers look.

## Importing from quantum chemistry

`g2spinach` is the central converter from a parsed electronic structure log to
`sys`/`inter`:

```matlab
[sys,inter]=g2spinach(props,particles,references,options)
```

- `props` — output of `gparse`, `oparse` or a compatible parser. EPR import requires explicit HFC source isotopes in `props.isotopes` (Gaussian mass numbers or ORCA isotope strings); missing or malformed provenance for selected, nonempty tensors is rejected before processing or warnings. Whole tensors are scaled by the target/source gyromagnetic-ratio ratio before thresholding and purging, while NMR import is unchanged. Empty tensors and electron-only selections need no provenance; zero-gamma sources are rejected, while direct zero-spin targets produce zero tensors.
- When selecting nuclei from a parsed EPR log, apply the same atom indices
  to `props.symbols`, `props.hfc.full.matrix`, `props.isotopes`, and
  `props.std_geom` when coordinates are included. If `std_geom` is absent,
  select from `props.inp_geom`; alternatively set `options.no_xyz=1` to omit
  coordinates. For a selected nonempty HFC, an isotope-array length mismatch
  now gives an alignment diagnostic; an unsupported source isotope identifies
  the atom and asks you to check isotope/symbol order. A known isotope
  without tabulated spin data gets a distinct atom-specific diagnostic. Empty
  or unselected tensors and electron-only or NMR imports need no source-isotope
  alignment.

- `particles` — cell array of element/isotope pairs, e.g.
  `{{'C','13C'},{'N','15N'}}`. Including an electron, as in
  `{{'E','E'},{'H','1H'}}`, switches the function into **EPR mode**: chemical
  shielding and scalar couplings are ignored, g-tensor and hyperfine couplings
  are imported, and coordinates are not returned.
- `references` — absolute isotropic shieldings of the reference substances, one
  per entry in `particles`, computed at the same level of theory. The header of
  `g2spinach.m` tabulates TMS absolute shieldings for GIAO and CSGT with B3LYP
  and HF at 6-31G* and 6-311+G(2d,p). Ignored in EPR mode.
- `options.min_j` — J-coupling threshold in Hz; `options.min_hfc` — hyperfine
  threshold in Hz on the Frobenius norm; `options.purge='on'` removes spins
  below `min_hfc` in EPR mode; `options.no_xyz=1` keeps only the interaction
  tensors and discards coordinates.

Outputs are `sys.isotopes`, `inter.coordinates` (cell array, angstrom),
`inter.zeeman.matrix` (ppm for nuclei, g-tensor for electrons),
`inter.coupling.matrix` and `inter.coupling.scalar` (both Hz), and
`inter.spinrot.matrix`.

```matlab
% NMR mode
[sys,inter]=g2spinach(gparse('../standard_systems/glycine.log'),...
                    {{'C','13C'},{'N','15N'}},[182.1 264.5],[]);
sys.magnet=14.1;

% EPR mode
options.no_xyz=1;
[sys,inter]=g2spinach(gparse('../standard_systems/chrysene_cation.log'),...
                            {{'E','E'},{'H','1H'}},[0 0],options);
sys.magnet=3.5;
```

`gparse(filename,options)` reads Gaussian 03/09/16 logs and returns geometries
in angstrom, energies in hartree, hyperfine tensors in gauss, J- and
K-couplings in Hz, plus `g`, `cst` (absolute shielding), `srt`, `nqi` and
`chi`. Options `'g_nosymm'`, `'cst_nosymm'`, `'hfc_nosymm'` disable tensor
symmetrisation.

`oparse(file_name)` reads ORCA 2.6–6.1 text output: `std_geom` in angstrom,
energy in hartree, ZFS in cm⁻¹, `hfc` in gauss, `efg` in a.u.⁻³, `nqi` in Hz,
`cst` in ppm, J-couplings in Hz, `chi` in cm³·K/mol with `chi_temps` in kelvin.
Its `props` can be handed to `g2spinach` or mined directly, as in
`props=oparse('cu_porph_hfc.out'); hfcs=props.hfc.full.matrix(26:37);`.
`ocparse(filename,pad_factor)` reads ORCA spin-density cube files in "3D simple
format".

`c2spinach(file_name)` reads the `[atoms]` and `[magres]` blocks of a CCP-NC
magres v1.0 file (CASTEP, Quantum ESPRESSO GIPAW) and returns `std_geom`
(angstrom), `symbols`, `natoms`, and, when the file has them, `cst` (shielding
relative to the bare nucleus in vacuum, ppm, in the printed component order),
`efg` (a.u.), and `k_couplings` (isotropic reduced couplings from the `isc`
records, in the same units as `gparse`, so that `g2spinach` converts them into
J-couplings for the isotopes it is given). Tensors are matched to atoms by
label and index, so a file whose `ms` records are reordered or partial still
lands on the right atoms; atoms without a tensor get an empty cell. Test for
optional fields with `isfield`. CASTEP shieldings must be referenced by hand,
and EFGs converted:

```matlab
props=c2spinach('mhc.magres');
inter.zeeman.matrix{n}=29.25*eye(3)-props.cst{n};
inter.coordinates={props.std_geom(2,:);
                   props.std_geom(5,:)};
nqi=castep2nqi(props.efg{5},20.44e-3,1);
inter.coupling.matrix{2,2}=remtrace(nqi);
```

`shift_iso(tensors,spin_numbers,new_iso)` replaces the isotropic parts of
imported shielding tensors with experimental isotropic shifts while keeping the
computed anisotropies, which is the standard way to combine DFT anisotropy with
measured shifts: `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,[1 2],
[43.6 110.0]);`.

`brokensymm(props_sing,props_trip)` estimates an exchange coupling in Hz from a
broken-symmetry singlet/triplet pair of DFT logs using the Yamaguchi equation
with the convention H = -2J(Sa·Sb). `karplus_fit(dir_path,atoms)` fits a
Karplus curve to a Gaussian dihedral scan.

## Importing structures and databases

`[sys,inter,aux]=protein(pdb_file,bmrb_file,options)` builds a protein spin
system from PDB coordinates and BMRB chemical shifts, guessing J-couplings from
Karplus curves and literature values and CSAs from local geometry.
`options.select` is `'backbone'`, `'backbone-minimal'`, `'backbone-hsqc'`,
`'all'`, or a list of PDB atom numbers; `options.pdb_mol` selects a molecule
from a multi-molecule PDB; `options.noshift` is `'keep'` (unassigned atoms
placed between -1 and 0 ppm) or `'delete'`; `options.deuterate` is a cell array
of PDB identifiers, or `'non-Me'`; `options.nh_csa` selects the peptide bond
CSA set, `'bax'`, `'tcb'` (default without a CSA file) or `'pol'`.
`options.csa_file` imports canonical AFNMR traceless tensors in ppm (five
comment lines, then atom header `serial name resname resnum` and three
matrix rows per atom). Every retained atom must match its PDB serial and
labels. Full Gaussian-row matrices are preserved, and BMRB scalar shifts
are unchanged. With a file, only an explicitly supplied `nh_csa` triggers
guessing: it warns and overwrites available amide N/H tensors, keeping
other imported anisotropies. Outputs include `sys.labels`
with IUPAC atom labels, `inter.coordinates` in angstrom, `inter.zeeman.scalar`
and `inter.zeeman.matrix` in ppm, `inter.coupling.scalar` in Hz, and `aux` with
residue numbers and types.

`[sys,inter]=nuclacid(pdb_file,shift_file,options)` does the same for nucleic
acids, taking shifts from an ASCII file formatted as
`[residue_number atom_id shift]`. `options.deut_list` marks deuterated atoms and
reduces the affected J-couplings; `options.noshift` behaves as above.
`read_pdb_pro(pdb_file_name,mod_id)`, `read_pdb_nuc(pdb_file_name)` and
`read_bmrb(bmrb_file_name)` are the underlying record readers.

`[sys,inter]=gissmo2spinach(filename,subsystem)` reads a GISSMO XML file and
returns a ready-to-use liquid-state NMR spin system. GISSMO supplies only
chemical shifts, J-couplings, a non-selective line width and the magnet field;
everything else has to be added by hand. The linewidth is Lorentzian FWHM
in hertz, converted to `pi*FWHM` inverse seconds. Its pure damping uses
`inter.rlx_keep='labframe'`: damping is added after retention, preserving
the spherical-tensor generator and supporting `zeeman-liouv` without
requesting unsupported diagonal retention.
`[sys,inter]=x2spinach(filename,shielding_refs)` reads SpinXML files.

Three importers handle data rather than parameters, and none of them produces
`sys` or `inter`. `vdata=v2spinach(inpath)` reads experimental FIDs in Varian
format from a data directory, returning `vdata.fid`, the acquisition header
fields, the fully parsed `procpar`, and derived spectral and diffusion
parameters; it has nothing to do with VASP. `bdata=b2spinach(inpath)` does the
same for Bruker experiment directories: fid or ser data, the acquisition and
processing parameter files, the digital filter group delay, and the gradient
and delay lists when present. `mesh=comsol_import(comsol)` imports a COMSOL 2D
mesh for the `meshflow` context from `comsol.mesh_file` and `comsol.velo_file`,
with `comsol.crop` and `comsol.inactivate` controlling the retained region.
`conc_plot(spin_system,conc,obs)` then draws concentrations on that mesh as
vertical bars after `mesh_plot` has drawn it, colouring each cell by one
(phase), two (phase and amplitude), or three (phase, amplitude, and
longitudinal) observables; a zero peak amplitude gives zero saturation and a
constant longitudinal observable gives full value, so such inputs render
instead of producing NaN face colours.

Complete ready-made systems live in `etc/molecules/` (`strychnine(spins)`,
`cyprinol()`, `lactate(spins)`, `allyl_pyruvate(spins)`,
`fatty_acid(nprotons)`, `dac_reaction()`) and `etc/diamond_defects/` (fifteen
`[sys,inter]=<name>(parameters)` builders for NV, P1 and related centres).
`guess_j_pro`, `guess_j_nuc` and `guess_csa_pro` are the estimators `protein`
and `nuclacid` call internally; everything they return is an estimate and must
be reported as one.

## Direct-sum operator addressing

In `sphten-liouv`, a product operator acts only in the substance hosting its
spins; cross-substance product specifications raise
`Spinach:which_subst:crossSubstance`. Isotope selections sum single-spin
operators across the hosting blocks. Operator and identity dimensions come
from `bas.offsets(end)`, not from the number of descriptor cells.
Explicit identity requests act only on the selected substance; numeric and
isotope sums contribute once per matching spin. Left/right identity actions
give the local identity, anticommutators twice it, and commutators zero.
Multi-substance Zeeman operators are local tensor products placed into their hosting direct-sum blocks. Isotope sums cover all matching substances.

Symmetry factorisations live in `bas.sym_fact(n)`. `reduce` and `rspt_eig`
embed their local projector columns using `bas.offsets`; do not read the
retired global `bas.irrep` field.

Converters preserve the substance direct sum. `sphten2zeeman` maps each
unit coordinate to `vec(I_D)` (Hilbert trace divided by `D` is the source
unit coordinate). `hilb2liouv` accepts explicit cells of Hilbert blocks;
numeric matrices retain their single-block meaning. `sim2liouv` uses
compiled block dimensions, migrates `sym_fact`, and refreshes the hash
from the descriptor cells, new dimensions, and substance membership.

Spatial phantom operators and states use the entire substance direct sum.
`imaging` and `meshflow` obtain its dimension from `bas.offsets(end)`;
`gridfree` checks extra isotropic terms against the same dimension.

Spin-only symmetry projectors are not applied to enlarged spatial-spin
generators in `reduce`; those use the usual trajectory-level reductions.
`v2fplanck` tensors the spatial transport with the full direct-sum identity.

`summary_basis` reports substances separately; `stateinfo` identifies the
hosting substance and global spin labels of every reported coefficient.
State numbers retain the global direct-sum offsets.

`zte` preserves all compiled spherical-tensor unit coordinates, even for
zero-population substances. `reduce` transports their support through the
symmetry projection and keeps it during subsequent population screening.

Trajectory analysis uses local descriptors and global spin labels. `trajan`
removes each unit independently and converts level populations per substance;
`trajsimil` groups equivalent tracks only within the same substance.

`kill_spin` and `dilute` rebuild an existing basis after particle removal.
Local manual columns and global numeric filter labels are reindexed. Isotope
filters are removed when no matching spin survives in their own substance;
empty substances keep their unit coordinate. Retained `inter_level`, `prox_level`,
and `space_level` depths are capped by each substance's surviving particle count
before the rebuild. IK-DNP vector components also respect electron and nucleus
counts; IK-SBS components respect mode, total active-particle, and spin counts.
Losing a required class in IK-DNP or IK-SBS selects IK-0 at the surviving
class depth (at least one), without connectivity pruning. Symmetry and Hamiltonian assumptions
are cleared, so reapply `assume` before constructing a Hamiltonian.

Cross-substance pair couplings raise `Spinach:create:crossSubstanceCoupling`;
product operators and states raise `Spinach:which_subst:crossSubstance`.

Chemistry uses explicit `chem.reactions` records through `kinetics` and
`react_gen`; retired rates, flux, and radical-pair fields are rejected.

Hilbert `evolution` reads the per-substance approximation cell; every Hilbert
block must use `none`, as enforced by `basis`.

`sim2liouv` accepts sparse horizontal density-matrix stacks and preserves their
column order while extracting each substance block. Segmented generators,
operator-like parameters, and state stacks must be block diagonal in substance;
nonzero cross-substance entries raise `Spinach:sim2liouv:crossSubstance`.

Synthetic compiled-system fixtures must supply offsets and local descriptor
cells too; bypassing `basis` does not restore the retired global layout.

`bootstrap` follows the same one-cell approximation contract as physical systems.

Two-substance chemistry fixtures require two approximation cells. The
generator and invariant suites cover positive matched multi-substance
reaction maps, explicit first-order exchange, routing, and conservation,
alongside single-substance spin permutations and empty reaction maps.

Single-substance Zeeman symmetry remains available through `bas.sym_fact(1)`;
multi-substance Zeeman symmetry and analytical filters remain explicitly rejected.
Zeeman units and equilibrium states are assembled independently per substance;
geometric units retain stock normalisations, while thermal blocks have trace c_n.
Single-substance descriptor consumers use `bas.basis{1}` and dimension consumers
use `bas.offsets(end)`. Identity states retain the selected substance: each
selected spin contributes one local unit in a sum, a local product contributes
once, and every `state` method weights that unit by its hosting concentration.
Use `coil_state` with the same description for unweighted detection vectors;
`state(...,'chem')` is a deprecated alias of the weighted exact method.

Imaging tests the compiled symmetry projectors, so a declared group disabled
through `sys.disable` does not prevent an imaging calculation.

Per-substance `space_level` aliases are derived independently; an empty entry
does not inherit the preceding substance's proximity depth.

Segmented Zeeman-Liouville states use local identities. The weighted `state`
wrapper rejects segmented wavefunctions with `Spinach:state:segmentedZeeman`;
`coil_state` can store the unweighted direct sum of local product kets.
Single-substance state numerics are unchanged.

Before basis compilation, `chem.parts` must cover every global spin; omitted
spins raise `Spinach:basis:incompletePartition`. Empty substances are permitted.

Level projectors retain their substance through mixed identity and non-identity
tensor expansions, including chemical concentration weighting.

`partner_state` keeps global descriptor positions but constructs states only in
the substance hosting the fixed and partner spins. Padding identities in its
returned descriptors do not imply populations in other substances.

Segmented coherent states and Zeeman steady solves are explicitly deferred.
Steady-state `solid_effect` and both DNP scans require a single substance;
these experiments retain their supported single-substance algorithms.

### Explicit chemistry records

Use `inter.chem.reactions`, a cell array of records containing `reactants`, `products`, `matching`, and `rate`. Substance indices are row vectors (repeats carry stoichiometry); matching is a two-column global spin map, and `zeros(0,2)` is an empty map. Empty products denote untracked loss. Rates may be non-negative scalars or time handles. `closure` defaults to `additive`; select `product` explicitly to retain cross-reactant polarisation products. Legacy rates/flux/radical-pair input fields are retired. Named selectors carry two electron indices on a single reactant; user selector matrices are substance-local. `merge_inp` shifts record substance and spin indices, not local selector matrices; entirely empty chemistry groups retain the default single-substance contract.

`unit_state`, `state`, and `equilibrium` return concentration-weighted density
matrices and Liouville states; storage-only wavefunctions remain unweighted.
Geometric detection and normalised operator vectors use `coil_state`, including
pulse-sequence receivers, voxel projection operators, normalised ENDOR sums,
relaxation-analysis probes, and microwave transition operators. Prepared densities
retain `state`; do not apply concentration factors to receivers. IME
`thermalize` instead takes unit-concentration target shapes: request equilibrium
on a copy with all `chem.concs` entries one, as `relaxation` does internally.
The propagated unit coordinates supply the instantaneous concentrations; neither
thermalisation nor pumping divides by a concentration. `magpump` takes an
unweighted `coil_state` target, and `steady` pins the supplied concentrations.

### Reaction propagation

`kinetics` compiles reaction records once. Numeric first-order records give a sparse matrix (zero numeric higher-order records do not change that classification); mass action and time-rate records give `K(t,eta)`. Time-only rate callbacks are evaluated once per assembly and shared across voxels. Use `1i*K` in a Liouvillian, or the existing `step` handle route for nonlinear propagation. `chem_concs` reads per-voxel concentrations from spherical-tensor unit coordinates or Zeeman trace functionals without division or normalisation. Spin-free pools participate dynamically. `react_gen` returns product-row/source-index lists, not the retired per-reactant generator matrices. Matched repeated spin-bearing reactants or products require occurrence-resolved matching and are rejected rather than assigned arbitrary molecular copies.

Use an explicit scalar `nz_shift`; the old `'chem'` shorthand is not defined for a general reaction network. `kill_spin` rebuilds reaction matching and basis data, but refuses removal of selector electrons or changes to a substance carrying user-supplied selector matrices.

### Intermolecular spin replacement

For a molecule A exchanging one spin with a pool B, use an additive `A+B -> A+B` record with matching that swaps those spins and retains the others. Departing-spin intramolecular correlations are destroyed; unaffected internal orders are retained. The concentrations are invariant because both sides have identical stoichiometry. With time-independent rates and no other concentration-changing reactions, evaluate the returned handle once at `unit_state(spin_system)` to obtain the constant additive generator for ordinary linear propagation. Do not freeze general mass-action or product-closure chemistry this way. `relayed_hyperpol` applies this construction to one ten-proton peptide block and twenty independent water pools; only water spins 11–20 participate in its forty replacement records.

`uf6_collisions` uses two first-order records with both nuclei matched in each direction; stationary concentrations weight excitation once, and `coil_state` detects both species without another population factor.

Worked replacements are `flux_asymmetric`/`flux_symmetric` and `frydman_pump_a`/`frydman_pump_b`: store each independent pool as its own substance, use one symmetric replacement per distinct pair, and use unweighted detection vectors. A one-spin pool has a complete level-one basis. Changes in thermal preparation must be checked separately from reaction-map equivalence.

`plain_reaction` demonstrates an entirely spin-free network: construct a ghost seed with explicit substance records, then trace it with `kill_spin` before calling `kinetics`. All five remaining coordinates are concentrations.

Small kernel-path demonstrations are `bimolecular_closures`, `spinless_sink_network`, and `cidnp_transport`; their corresponding registered tests cover mass action, the two closures, selective loss, and integrated nuclear product arrival.

Reaction-bearing systems bypass spin-only symmetry factorisation in `reduce`: chemical maps can connect substance irreps. Full-generator ZTE and path tracing remain available and retain chemical arrival into initially empty products.

### Two-stage chemistry histories

`diels_alder_zmag`, `diels_alder_spec`, and `reacting_nmr` put true initial concentrations in `chem.concs`, trace spins for their LG4 concentration histories, and compile the full additive maps once for the two-point spin steps. Embed the prescribed history into unit coordinates when evaluating `K(t,eta)`; initialise spin magnetisation with weighted `state`, add `unit_state`, and detect with unweighted `coil_state`. Do not multiply the initial state by the concentrations a second time.

### Spatial two-stage chemistry

`reacting_flow_nmr` traces all spins with `kill_spin` for its concentration-only stage: the resulting five unit blocks include the traced solvent, and the same reaction records generate both concentration and spin transport. Its frozen-rate stepping workflow and `makima` history are retained, but equal sharing of product unit arrival changes finite frozen concentration steps relative to the old asymmetric generator. Equal sharing is the specified additive closure, not a claim of improved frozen-step accuracy: both allocations converge to the same mass-action ODE. Do not claim numerical history equivalence from the equal instantaneous derivative. In the NMR stage, history values are placed into voxel unit coordinates before evaluating `K(t,eta)`; additive closure then depends only on those coordinates, not the spin orders. The actual propagated state includes unit populations, while detection and reference longitudinal vectors use `coil_state`. The shared `dac_reaction` definition retains all three acetonitrile protons and their T1/T2 relaxation. The flow example prepares only the four reacting species; solvent remains unexcited and contributes no signal.

`create(sys)` without interaction input remains supported; it uses an empty
interaction structure and the ordinary one-substance/unit-concentration defaults.

In chemistry-free systems, `reduce` rejects cross-substance entries in caller-supplied generators at the
compiled spin dimension before building substance-local projectors. This
boundary also covers the adjoint generator passed by destination screening.

The unweighted primitive requires `coil_state(spin_system,states,spins,method)`
with all four arguments; use `exact` or `cheap`, and pass `[]` for wavefunction
spin lists. Only the legacy `state` wrapper retains optional arguments.

Zeeman-wavefunction direct sums provide storage only. Nonempty reaction records
are rejected by `basis` and `kinetics`; `equilibrium`, `unit_state`, and
`thermalize` reject mixed-state or concentration-weighted requests explicitly
with messages naming `zeeman-wavef`. Legacy unweighted singleton kets remain available.

### Zeeman chemistry and matrix IME

Zeeman Liouville chemistry uses complete local bases and dense transformations
of the shared reaction compiler, so it is intended as a small-system reference
backend. Both closures, atom matching, spin-free products, selectors, and spatial
state-dependent maps retain the same physical contract. User selector pairs
remain local Liouville product superoperators even for Hilbert input.

With Hilbert reaction records, `K=kinetics(s)` returns `K(t,rho)` as a matrix
derivative, not a conjugation Hamiltonian. Only first-order records are accepted;
mass action raises `Spinach:kinetics:hilbertMassAction`. `thermalize` supports
Hilbert IME from an explicit direct-sum Liouville relaxation matrix and a
block-diagonal unit-trace target, returning the same matrix-RHS interface. These
handles are not `step` generators, and the stock Hilbert relaxation constructor
is unchanged. `chem_concs` extracts each Hilbert block trace and rejects
inter-substance coherences.

`unit_state` requires a compiled basis; absent basis metadata is rejected before formalism capability checks.

Retired global basis matrices and `bas.irrep` are rejected at the `basis`,
`coherence`, `correlation`, `summary_basis`, and `kinetics` boundaries. Use
`bas.basis{n}` with `bas.offsets`, and `bas.sym_fact(n).irr_projectors`/
`irr_dimensions`. `kinetics` also rejects retired chemistry fields inserted
after `create`; use `chem.reactions`. These are consumer checks, not custom
dot-read interception: the compiled objects remain ordinary MATLAB structs.

Linear contexts (`liquid`, `imaging`, `crystal`, `powder`, `device`, `floquet`, `singlerot`, `doublerot`, and `gridfree`) reject state-dependent kinetics handles. Use a custom pulse sequence with `step`/`iserstep` for multi-reactant or callback-rate records; `examples/kinetics/nonlinear/bimolecular_closures.m` and `examples/microfluidics/reacting_flow_nmr.m` demonstrate that route.
