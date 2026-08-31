# experiments/nmr_liquids/gcosy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/gcosy.m`
- Signature: `fid=gcosy(spin_system,parameters,H,R,K)`
- Total lines: 266

## Purpose

Horne-Morris gradient-selected COSY pulse sequence. Syntax: fid=gcosy(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 59-60: Set default gradient amplitude; implemented by `if ~isfield(parameters,'g_amp')`.
- Lines 64-65: Set default gradient duration; implemented by `if ~isfield(parameters,'g_dur')`.
- Lines 69-70: Set default gradient stabilisation delay; implemented by `if ~isfield(parameters,'g_stab_del')`.
- Lines 74-75: Set default active sample length; implemented by `if ~isfield(parameters,'s_len')`.
- Lines 79-80: Set default coherence pathway; implemented by `if ~isfield(parameters,'pathway')`.
- Lines 84-85: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 87-88: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 90-91: Compute the evolution timestep; implemented by `timestep=1/parameters.sweep`.
- Lines 93-94: Initial state up to a constant multiplier; implemented by `rho=state(spin_system,'Lz',parameters.spins{1},'cheap')`.
- Lines 96-97: Detection state up to a constant multiplier; implemented by `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 99-100: Get the pulse operator; implemented by `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 102-103: Apply the first pulse; implemented by `rho=step(spin_system,Lx,rho,pi/2)`.
- Lines 105-106: Run the F1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,timestep,parameters.npoints(1)-1,'trajectory')`.
- Lines 108-109: Build the second pulse propagator; implemented by `P=propagator(spin_system,Lx,parameters.angle)`.
- Lines 111-112: Include the recovery delay after the first gradient; implemented by `if parameters.g_stab_del>0`.
- Lines 116-117: Select the requested gradient pathway; implemented by `switch upper(parameters.pathway)`.
- Lines 121-123: Use opposite signs for P-type selection; implemented by `rho_stack=grad_sandw(spin_system,L,rho_stack,P,parameters.g_amp*[1 -1], parameters.s_len,parameters.g_dur*[1 1],[1 1])`.
- Lines 125-126: Run the recovery delay after the second gradient; implemented by `if parameters.g_stab_del>0`.

### Control flow inferred from the code

- Line 60: conditional branch on `~isfield(parameters,'g_amp')`.
- Line 65: conditional branch on `~isfield(parameters,'g_dur')`.
- Line 70: conditional branch on `~isfield(parameters,'g_stab_del')`.
- Line 75: conditional branch on `~isfield(parameters,'s_len')`.
- Line 80: conditional branch on `~isfield(parameters,'pathway')`.
- Line 112: conditional branch on `parameters.g_stab_del>0`.
- Line 117: dispatches on `upper(parameters.pathway)`; cases `'P'`, `'N'`, `'P+N'`.
- Line 126: conditional branch on `parameters.g_stab_del>0`.
- Line 141: conditional branch on `parameters.g_stab_del>0`.
- Line 160: conditional branch on `parameters.g_stab_del>0`.

### Key state/data transformations

- Lines 61: computes `parameters.g_amp` using `parameters.g_amp=3`.
- Lines 66: computes `parameters.g_dur` using `parameters.g_dur=2e-3`.
- Lines 71: computes `parameters.g_stab_del` using `parameters.g_stab_del=2e-4`.
- Lines 76: computes `parameters.s_len` using `parameters.s_len=1.5`.
- Lines 81: computes `parameters.pathway` using `parameters.pathway='P'`.
- Lines 88: computes `L` using `L=H+1i*R+1i*K`.
- Lines 91: computes `timestep` using `timestep=1/parameters.sweep`.
- Lines 94: computes `rho` using `rho=state(spin_system,'Lz',parameters.spins{1},'cheap')`.
- Lines 97: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 100: computes `Lx` using `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 106: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,timestep,parameters.npoints(1)-1,'trajectory')`.
- Lines 109: computes `P` using `P=propagator(spin_system,Lx,parameters.angle)`.
- Lines 131-132: computes `fid` using `fid=evolution(spin_system,L,coil,rho_stack,timestep, parameters.npoints(2)-1,'observable')`.
- Lines 152-153: computes `rho_stack_p` using `rho_stack_p=grad_sandw(spin_system,L,rho_stack,P,parameters.g_amp*[1 -1], parameters.s_len,parameters.g_dur*[1 1],[1 1])`.
- Lines 156-157: computes `rho_stack_n` using `rho_stack_n=grad_sandw(spin_system,L,rho_stack,P,parameters.g_amp*[1 1], parameters.s_len,parameters.g_dur*[1 1],[1 1])`.
- Lines 166-167: computes `fid.pos` using `fid.pos=evolution(spin_system,L,coil,rho_stack_p,timestep, parameters.npoints(2)-1,'observable')`.
- Lines 168-169: computes `fid.neg` using `fid.neg=evolution(spin_system,L,coil,rho_stack_n,timestep, parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 176: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Parameters / inputs

- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.angle second pulse angle in radians, usu-
- ally pi/2, but also allows COSY45,
- COSY60, etc.
- parameters.g_amp gradient amplitude in Gauss/cm,
- defaults to 3
- parameters.g_dur gradient duration in seconds,
- defaults to 2e-3
- parameters.g_stab_del post-gradient stabilisation delay in
- seconds, defaults to 2e-4
- parameters.s_len active sample length in cm,
- defaults to 1.5
- parameters.pathway optional coherence pathway selection,
- either 'P', 'N', or 'P+N', defaults
- to 'P'
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -two-dimensional free induction decay, or a structure
- with P-type fid.pos and N-type fid.neg fields in
- 'P+N' mode
- Note: the default P-type pathway uses opposite gradient signs
- and is less sensitive to mixing pulse phase errors. The
- N-type pathway uses equal gradient signs.
- Note: 'P+N' mode returns P-type and N-type components for
- echo/anti-echo recombination in phase-sensitive processing.

## Implementation structure

- Horne-Morris gradient-selected COSY pulse sequence. Syntax:
- fid=gcosy(spin_system,parameters,H,R,K)
- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.angle second pulse angle in radians, usu-
- ally pi/2, but also allows COSY45,
- COSY60, etc.
- parameters.g_amp gradient amplitude in Gauss/cm,
- defaults to 3
- parameters.g_dur gradient duration in seconds,

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isfield()`, `grumble()`, `state()`, `operator()`, `step()`, `evolution()`, `propagator()`, `upper()`, `grad_sandw()`, `ismember()`, `ismatrix()`, `all()`, `elseif()`, `iscell()`, `ischar()`, `any()`.
