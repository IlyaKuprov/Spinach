# kernel/eigenfields/eigenfields.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/eigenfields/eigenfields.m`
- Signature: `tran=eigenfields(spin_system,parameters,Hz,Hc,Hmw)`
- Total lines: 514

## Purpose

Computes resonance fields. For a Hamiltonian Hc+B*Hz, returns all magnetic fields B for which the difference between two eigenvalues of Hc+B*Hz is equal to the frequency provided, and the transition moment across the specified operator Hmw is significant. Syntax: tran=eigenfields(spin_system,parameters,Hz,Hc,Hmw)

## Physical / mathematical content

- Eigenfield utilities. These files analyse field-dependent eigenstructure and resonance conditions, linking Hamiltonian spectra to magnetic-field sweeps and transition behaviour.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 65-66: Check consistency; implemented by `grumble(spin_system,parameters,Hz,Hc)`.
- Lines 68-69: Get frequency gap tolerance in rad/s; implemented by `frq_gap_tol=abs(spin('E')*parameters.pp_tol)`.
- Lines 71-72: Cubic spline root tolerance; implemented by `root_tol=sqrt(eps)`.
- Lines 74-75: Get MW frequency in rad/s; implemented by `omega=2*pi*parameters.mw_freq`.
- Lines 77-78: Pathway selection; implemented by `switch spin_system.bas.formalism`.
- Lines 80-81: Hilbert space; implemented by `case 'zeeman-hilb'`.
- Lines 83-84: Initial field grid: trisect the window; implemented by `left_edge=min(parameters.window)`.
- Lines 91-92: Set up field grid data structures; implemented by `E= zeros([size(Hz,1) numel(grid)])`.
- Lines 98-99: Populate field grid; implemented by `for n=1:numel(grid)`.
- Lines 101-102: Sorted (ascending) energies, dE/dB derivatives, transition moments, level pops; implemented by `[E(:,n),V{n},dE(:,n),T{n},LP(:,n)]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,grid(n))`.
- Lines 106-107: Paint all intervals dirty; implemented by `converged=false(numel(grid)-1,1)`.
- Lines 109-110: Iterative trisection; implemented by `while any(~converged)`.
- Lines 112-113: Start at the leftmost point of the field grid; implemented by `new_grid=grid(1); new_E=E(:,1); new_dE=dE(:,1)`.
- Lines 116-117: Inspect old intervals; implemented by `for n=2:numel(grid)`.
- Lines 119-120: Check for prior convergence; implemented by `if converged(n-1)`.
- Lines 122-123: Inherit the flag; implemented by `new_conv(end+1)=true()`.
- Lines 127-128: Field grid interval; implemented by `dx=grid(n)-grid(n-1)`.
- Lines 130-133: Predict and sort midpoint energies; implemented by `EmP1=herm_spline(E(:,n-1),dE(:,n-1)*dx, E(:,n), dE(:,n)*dx, (1/3)*ones(size(E,1),1))`.

### Control flow inferred from the code

- Line 78: dispatches on `spin_system.bas.formalism`; cases `'zeeman-hilb'`.
- Line 99: `for` loop over `n=1:numel(grid)`.
- Line 110: `while` loop over `any(~converged)`.
- Line 117: `for` loop over `n=2:numel(grid)`.
- Line 120: conditional branch on `converged(n-1)`.
- Line 153: conditional branch on `(norm(EmP1-Em1,1)<frq_gap_tol)&&`.
- Line 183: conditional branch on `numel(grid)>1e3, error('field grid construction failed'); end`.
- Line 188: `for` loop over `n=2:numel(grid)`.
- Line 194: `for` loop over `k=1:size(E,1)`.
- Line 198: conditional branch on `max_ovlp<sqrt(eps), error('state tracking failed'); end`.
- Line 206: conditional branch on `any(~state_hit), error('state tracking failed'); end`.
- Line 217: `for` loop over `n=1:numel(T)`.
- Line 232: `for` loop over `n=2:numel(grid)`.
- Line 238: conditional branch on `isempty(source), continue; end`.

### Key state/data transformations

- Lines 69: computes `frq_gap_tol` using `frq_gap_tol=abs(spin('E')*parameters.pp_tol)`.
- Lines 72: computes `root_tol` using `root_tol=sqrt(eps)`.
- Lines 75: computes `omega` using `omega=2*pi*parameters.mw_freq`.
- Lines 84: computes `left_edge` using `left_edge=min(parameters.window)`.
- Lines 85: computes `right_edge` using `right_edge=max(parameters.window)`.
- Lines 86-89: computes `grid` using `grid=[left_edge+0*(right_edge-left_edge)/3 left_edge+1*(right_edge-left_edge)/3 left_edge+2*(right_edge-left_edge)/3 left_edge+3*(right_edge-left_edge)/3]`.
- Lines 92: computes `E` using `E= zeros([size(Hz,1) numel(grid)])`.
- Lines 93: computes `LP` using `LP=zeros([size(Hz,1) numel(grid)])`.
- Lines 94: computes `dE` using `dE=zeros([size(Hz,1) numel(grid)])`.
- Lines 95: computes `V` using `V=cell(1,numel(grid))`.
- Lines 96: computes `T` using `T=cell(1,numel(grid))`.
- Lines 102: computes `[E(:,n),V{n},dE(:,n),T{n},LP(:,n)]` using `[E(:,n),V{n},dE(:,n),T{n},LP(:,n)]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,grid(n))`.
- Lines 107: computes `converged` using `converged=false(numel(grid)-1,1)`.
- Lines 113: computes `new_grid` using `new_grid=grid(1); new_E=E(:,1); new_dE=dE(:,1)`.
- Lines 114: computes `new_T` using `new_T=T(1); new_V=V(1); new_LP=LP(:,1); new_conv=[]`.
- Lines 123: computes `new_conv(end+1)` using `new_conv(end+1)=true()`.
- Lines 128: computes `dx` using `dx=grid(n)-grid(n-1)`.
- Lines 131-133: computes `EmP1` using `EmP1=herm_spline(E(:,n-1),dE(:,n-1)*dx, E(:,n), dE(:,n)*dx, (1/3)*ones(size(E,1),1))`.

### Local helper functions

- Line 444: `grumble()` — `function grumble(spin_system,parameters,Hz,Hc)`.
  - Representative operation: `if (~isnumeric(Hc))||(size(Hc,1)~=size(Hc,2))|| (~isnumeric(Hz))||(size(Hz,1)~=size(Hz,2))|| (size(Hc,1)~=size(Hz,1))`.
  - Representative operation: `(~isnumeric(Hz))||(size(Hz,1)~=size(Hz,2))|| (size(Hc,1)~=size(Hz,1))`.

## Parameters / inputs

- Hz -field-dependent part of the laboratory frame Hamil-
- tonian operator (Hilbert space) or commutation supe-
- roperator (Liouville space), normalised to 1 Tesla
- Hc -field-independent part of the laboratory frame Hamil-
- tonian operator (Hilbert space) or commutation supe-
- roperator (Liouville space), containing couplings and
- offsets
- Hmw -observable operator (Hilbert space) or observable
- vector (Liouville space), without the amplitude
- prefactor
- parameters.window -magnet field window, Tesla
- parameters.mw_freq -microwave frequency, Hz
- parameters.orientation -three Euler angles in radians
- specifying the system orientation
- parameters.tm_tol -relative transition moment
- tolerance
- parameters.pp_tol -peak position tolerance in Tesla,
- this should be much smaller than
- the typical line width
- parameters.fwhm -transition full width at half
- maximum, Tesla
- parameters.rspt_order -perturbation theory order to use
- to account for the off-diagonal
- part of the Hamiltonian, Inf for
- exact diagonalisation

## Outputs

- tran.tf -vector of transition fields in Tesla
- tran.tm -vector of transition moments
- tran.tw -vector of transition FWHMs in Tesla
- tran.pd -vector of energy level population differences
- tran.ti -transition identity array, one row per transition
- tran.tj -vector of scaled field-sweep Jacobians

## Implementation structure

- Computes resonance fields. For a Hamiltonian Hc+B*Hz, returns all
- magnetic fields B for which the difference between two eigenvalues
- of Hc+B*Hz is equal to the frequency provided, and the transition
- moment across the specified operator Hmw is significant. Syntax:
- tran=eigenfields(spin_system,parameters,Hz,Hc,Hmw)
- Hz - field-dependent part of the laboratory frame Hamil-
- tonian operator (Hilbert space) or commutation supe-
- roperator (Liouville space), normalised to 1 Tesla
- Hc - field-independent part of the laboratory frame Hamil-
- roperator (Liouville space), containing couplings and
- offsets
- Hmw - observable operator (Hilbert space) or observable

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `rspt_eig()`, `false()`, `any()`, `converged()`, `new_conv()`, `true()`, `herm_spline()`, `new_grid()`, `new_dE()`, `new_LP()`, `new_E()`, `new_T()`, `new_V()`, `state_ovlp()`.
