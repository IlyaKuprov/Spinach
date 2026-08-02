# Spinach input specification

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

Legal `sys.disable` entries, anything else being an error: `'zte'`
(unpopulated-state elimination), `'pt'` (non-interacting subspace detection),
`'symmetry'`, `'krylov'`, `'clean-up'`, `'hygiene'` (start-up health checks),
`'dss'`, `'expv'`, `'trajlevel'`, `'merge'`, `'colorbar'`, `'asyredf'`. Legal
`sys.enable` entries: `'gpu'`, `'op_cache'`, `'ham_cache'`, `'prop_cache'`,
`'xmemlist'`, `'greedy'`, `'paranoia'` (tight tolerances), `'cowboy'` (loose
tolerances), `'polyadic'`, `'sodd'` (spin-orbit corrections to dipolar
couplings), `'dafuq'`.

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
| `'E'` | Electron, multiplicity 2. |
| `'E4'`, `'E16'` | High-spin electron; the integer is the **multiplicity**, so `'E4'` is S=3/2 and `'E16'` is S=15/2. |
| `'G'` | Ghost spin, gamma=0 and multiplicity 1; a placeholder that carries coordinates but no magnetism. |
| `'N'`, `'M'` | Neutron, muon. |
| `'C#'`, `'V#'`, `'T#'` | Cavity mode, phonon mode, transmon; the integer is the number of levels, and gamma is zero for all three. |

`iselectron` and `isnucleus` test a specification string. `isoswap(sys,inter,
spins,new_iso)` performs isotope replacement and rescales all interactions
accordingly, wiping quadratic and higher-order couplings with a warning.

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
in `coeff` must have its Euler angles supplied.

```matlab
inter.giant.coeff={{[0 0 0],Bkq{2},[0 0 0 0 0 0 0],Bkq{4}}};
inter.giant.euler={{[0 0 0],[0 0 0],[0 0 0],[0 0 0]}};
```

Stevens operator coefficients are converted to the irreducible spherical tensor
convention with `stev2sph(k,Bkq)`.

## Chemistry and kinetics

```matlab
inter.chem.parts={1,2};                        % spin index sets per species
inter.chem.rates=[-2e4   2e4
                   2e4  -2e4];                 % Hz, columns sum to zero
inter.chem.concs=[1.0 1.0];
```

- `parts` — cell array of index vectors, one per chemical subsystem, disjoint
  and within the spin count. The default is one subsystem containing everything.
- `concs` — initial concentrations, one per subsystem, non-negative; mandatory
  as soon as there is more than one subsystem.
- `rates` — square first-order rate matrix in hertz, one row and column per
  subsystem, column sums negligible (conservation of matter is checked). In
  exchange mode all subsystems must have the same number of spins and identical
  isotope sequences, so that spin *k* of species A maps onto spin *k* of B.
- `flux_rate` and `flux_type` — magnetisation flux between individual spins;
  an `nspins x nspins` real matrix and either `'intermolecular'` or
  `'intramolecular'`. Both must be supplied together.
- `rp_theory`, `rp_electrons`, `rp_rates` — radical pair recombination:
  `'haberkorn'`, `'jones-hore'` or `'exponential'`; the two recombining
  electron indices; and `[singlet_rate triplet_rate]` in hertz. All three
  fields must appear together.

Any chemistry at all forces `bas.formalism='sphten-liouv'`.
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
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=3;
```

| Field | Legal values | Notes |
|---|---|---|
| `formalism` | `'sphten-liouv'`, `'zeeman-liouv'`, `'zeeman-hilb'`, `'zeeman-wavef'` | Mandatory. |
| `approximation` | `'none'`, `'IK-0'`, `'IK-1'`, `'IK-2'`, `'IK-DNP'` | Mandatory. Only `'none'` is legal outside `sphten-liouv`. |
| `connectivity` | `'scalar_couplings'`, `'full_tensors'` | Required by, and only legal for, `IK-1` and `IK-2`. |
| `level` | positive integer; `1x3` integer vector for `IK-DNP` | Required by `IK-0`, `IK-1`, `IK-DNP`. Cannot exceed the number of spins. For `IK-DNP` the three entries bound electrons, spins and nuclei respectively. |
| `space_level` | positive integer | Required by, and only legal for, `IK-1` and `IK-2`. |
| `projections` | row vector of integers | Keeps only the listed total projection quantum numbers. `sphten-liouv` only. |
| `longitudinals` | cell array of isotope strings or spin index vectors | Keeps only longitudinal states on those spins. |
| `zero_quantum` | cell array of isotope strings or spin index vectors | Keeps only zero-quantum coherences on those spins. |
| `manual` | logical matrix with `nspins` columns | Explicit basis state list, one state per row. |
| `sym_group` | cell array from `S2`, `S3`, `S4`, `S4A`, `S5`, `S6`, `S6A`, `S8A` | Permutation symmetry groups. |
| `sym_spins` | cell array of index vectors | One vector per group, at least two spins each, no spin in two groups. Mandatory alongside `sym_group`. |
| `sym_a1g_only` | logical | Keep only the fully symmetric irreducible representation. |

`sphten-liouv` refuses multiplicities above 16. `IK-DNP` requires both
electrons and nuclei and nothing else in the system.

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

- `props` — output of `gparse`, `oparse` or a compatible parser.
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

`c2spinach(file_name)` reads the new-format section of a CASTEP `.magres` file
and returns `std_geom` (angstrom), `symbols`, `cst` (shielding relative to the
bare nucleus in vacuum, ppm) and `efg` (a.u.⁻³). CASTEP shieldings must be
referenced by hand, and EFGs converted:

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
CSA set, `'bax'`, `'tcb'` (default) or `'pol'`. Outputs include `sys.labels`
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
everything else has to be added by hand.
`[sys,inter]=x2spinach(filename,shielding_refs)` reads SpinXML files.

Two importers handle data rather than parameters, and neither produces `sys` or
`inter`. `vdata=v2spinach(inpath)` reads an experimental FID in Varian format
from a data directory, returning `vdata.fid` and the acquisition header fields;
it has nothing to do with VASP. `mesh=comsol_import(comsol)` imports a COMSOL 2D
mesh for the `meshflow` context from `comsol.mesh_file` and `comsol.velo_file`,
with `comsol.crop` and `comsol.inactivate` controlling the retained region.

Complete ready-made systems live in `etc/molecules/` (`strychnine(spins)`,
`cyprinol()`, `lactate(spins)`, `allyl_pyruvate(spins)`,
`fatty_acid(nprotons)`, `dac_reaction()`) and `etc/diamond_defects/` (fourteen
`[sys,inter]=<name>(parameters)` builders for NV, P1 and related centres).
`guess_j_pro`, `guess_j_nuc` and `guess_csa_pro` are the estimators `protein`
and `nuclacid` call internally; everything they return is an estimate and must
be reported as one.
