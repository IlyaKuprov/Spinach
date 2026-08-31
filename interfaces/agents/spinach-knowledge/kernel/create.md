# kernel/create.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/create.m`
- Signature: `spin_system=create(sys,inter)`
- Total lines: 3192

## Purpose

The entry function of the Spinach kernel that creates the spin system object that the rest of the library requires to run. It checks and absorbs interaction specifications, and writes some useful diagnostics to the console. Syntax: spin_system=create(sys,inter)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 47-48: Refuse to run in script mode; implemented by `call_stack=dbstack`.
- Lines 54-55: On head node, close all open files; implemented by `if ~isworkernode, fclose('all'); end`.
- Lines 57-58: Locate the root and run sanity checks; implemented by `if isempty(which('existentials'))`.
- Lines 60-62: Tell the user to RTFM; implemented by `error('paths have not been correctly set - please follow installation\n%s', 'instructions (make sure you had included the subdirectories).')`.
- Lines 66-67: Set the root directory; implemented by `root_dir=which('existentials')`.
- Lines 72-73: Overrides; implemented by `autoexec`.
- Lines 75-76: Rare, but it can happen; implemented by `if nargin==1, inter=[]; end`.
- Lines 78-79: Validate input; implemented by `grumble(sys,inter)`.
- Lines 81-83: Run existential checks; implemented by `if (~isfield(sys,'disable'))|| (~ismember('hygiene',sys.disable))`.
- Lines 87-88: Decide output destination; implemented by `if isfield(sys,'output')`.
- Lines 90-91: String specifications; implemented by `if strcmp(sys.output,'hush')`.
- Lines 93-94: Hush the output; implemented by `spin_system.sys.output='hush'`.
- Lines 98-99: Print to the console; implemented by `spin_system.sys.output=1`.
- Lines 103-104: Print to a user-specified file; implemented by `spin_system.sys.output=fopen(sys.output,'a')`.
- Lines 108-109: Parse out; implemented by `sys=rmfield(sys,'output')`.
- Lines 118-119: Decide scratch destination; implemented by `if isfield(sys,'scratch')`.
- Lines 121-122: Scratch to a user-specified directory; implemented by `spin_system.sys.scratch=sys.scratch`.
- Lines 124-125: Parse out; implemented by `sys=rmfield(sys,'scratch')`.

### Control flow inferred from the code

- Line 49: conditional branch on `strcmp(call_stack(end).name,'create')&&(~isworkernode)`.
- Line 55: conditional branch on `~isworkernode, fclose('all'); end`.
- Line 58: conditional branch on `isempty(which('existentials'))`.
- Line 76: conditional branch on `nargin==1, inter=[]; end`.
- Line 82: conditional branch on `(~isfield(sys,'disable'))||`.
- Line 88: conditional branch on `isfield(sys,'output')`.
- Line 91: conditional branch on `strcmp(sys.output,'hush')`.
- Line 119: conditional branch on `isfield(sys,'scratch')`.
- Line 135: conditional branch on `~exist(spin_system.sys.scratch,'dir')`.
- Line 146: conditional branch on `isfield(sys,'disable')`.
- Line 162: conditional branch on `isfield(sys,'enable')`.
- Line 181: conditional branch on `~isempty(spin_system.sys.disable)`.
- Line 183: conditional branch on `ismember('hygiene',spin_system.sys.disable), report(spin_system,' > health checks at start-up'); end`.
- Line 184: conditional branch on `ismember('pt',spin_system.sys.disable), report(spin_system,' > detection of non-interacting subspaces'); end`.

### Key state/data transformations

- Lines 48: computes `call_stack` using `call_stack=dbstack`.
- Lines 67: computes `root_dir` using `root_dir=which('existentials')`.
- Lines 68: computes `spin_system.sys.root_dir` using `spin_system.sys.root_dir=root_dir(1:(end-32))`.
- Lines 94: computes `spin_system.sys.output` using `spin_system.sys.output='hush'`.
- Lines 109: computes `sys` using `sys=rmfield(sys,'output')`.
- Lines 122: computes `spin_system.sys.scratch` using `spin_system.sys.scratch=sys.scratch`.
- Lines 149: computes `spin_system.sys.disable` using `spin_system.sys.disable=sys.disable`.
- Lines 165: computes `spin_system.sys.enable` using `spin_system.sys.enable=sys.enable`.
- Lines 178: computes `[spin_system,sys]` using `[spin_system,sys]=tolerances(spin_system,sys)`.
- Lines 214: computes `spin_system.sys.job_id` using `spin_system.sys.job_id=md5_hash([posixtime(datetime('now')) feature('getpid')])`.
- Lines 218-219: computes `spin_system.sys.job_dir` using `spin_system.sys.job_dir=[spin_system.sys.scratch filesep 'spinach_job_' spin_system.sys.job_id]`.
- Lines 246: computes `nworkers` using `nworkers=feature('numcores')-1`.
- Lines 249: computes `sys.parallel` using `sys.parallel={'processes',max([1 nworkers])}`.
- Lines 257: computes `current_pool` using `current_pool=gcp('nocreate')`.
- Lines 263: computes `delete(current_pool); current_pool` using `delete(current_pool); current_pool=gcp('nocreate')`.
- Lines 270: computes `pool_storage` using `pool_storage=current_pool.Cluster.JobStorageLocation`.
- Lines 272: computes `pool_job_dirs` using `pool_job_dirs=dir([pool_storage filesep 'Job*'])`.
- Lines 276-277: computes `value_store_dir` using `value_store_dir=[pool_job_dirs(n).folder filesep pool_job_dirs(n).name filesep 'value_store']`.

### Local helper functions

- Line 1525: `grumble()` — `function grumble(sys,inter)`. Check the output switch
  - Representative operation: `if isfield(sys,'output')`.
  - Representative operation: `if ~ischar(sys.output)`.

## Parameters / inputs

- sys -spin system and instrument specification
- structure, see the spin system specifica-
- tion section of the online manual
- inter -interaction specification structure, see
- see the spin system specification section
- of the online manual

## Outputs

- spin_system -the primary object used by Spinach
- to store simulation information
- Note: inter.modes.carriers declares, for each bosonic mode, the labo-
- ratory frequency of the rotating frame in which inter.modes.frqs
- is specified; declared frequencies are then detunings and may be
- negative, and thermal occupations are computed from the physical
- frequency, meaning the sum of the carrier and the detuning.
- Note: inter.modes.t2_times is interpreted at the declared temperature:
- the pure dephasing rate is extracted as 1/T2-kappa*(1+2*nbar)/2,
- where kappa is the amplitude damping rate and nbar the thermal
- occupation at the physical mode frequency.
- Note: quadrature operators of bosonic modes follow the (a+a')/sqrt(2)
- normalisation everywhere, in inter.modes.longitudinal as well as
- in the modulation channels.

## Implementation structure

- The entry function of the Spinach kernel that creates the spin
- system object that the rest of the library requires to run. It
- checks and absorbs interaction specifications, and writes some
- useful diagnostics to the console. Syntax:
- spin_system=create(sys,inter)
- sys -spin system and instrument specification
- structure, see the spin system specifica-
- tion section of the online manual
- inter -interaction specification structure, see
- see the spin system specification section
- of the online manual
- spin_system -the primary object used by Spinach

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strcmp()`, `call_stack()`, `fclose()`, `which()`, `instructions()`, `root_dir()`, `grumble()`, `isfield()`, `ismember()`, `existentials()`, `fopen()`, `rmfield()`, `exist()`, `mkdir()`, `report()`, `banner()`.
