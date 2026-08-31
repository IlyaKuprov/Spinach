# kernel/optimcon/bss_ops.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/bss_ops.m`
- Signature: `resp_ops=bss_ops(spin_system,channels,carrier_frq)`
- Total lines: 125

## Purpose

Bloch-Siegert response operators for the optimal control module. For each control channel, returns the operator whose coefficient in every time slice of a GRAPE optimisation is the square of the physical con- trol amplitude on that channel. The operator collects the second-order Bloch-Siegert frequency shifts of every spin in the system: B=sum_n (gamma_n/gamma_c)^2*[1/(2*(omega_n+omega_c)) +(foreign isotopes only) 1/

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 43-44: Check consistency; implemented by `grumble(spin_system,channels,carrier_frq)`.
- Lines 46-47: Preallocate the answer; implemented by `resp_ops=cell(size(channels))`.
- Lines 49-50: Loop over control channels; implemented by `for c=1:numel(channels)`.
- Lines 52-53: Get the channel magnetogyric ratio; implemented by `gamma_chan=spin(channels{c})`.
- Lines 55-56: Preallocate the response operator; implemented by `resp_op=mprealloc(spin_system,1)`.
- Lines 58-59: Loop over the spins; implemented by `for n=1:spin_system.comp.nspins`.
- Lines 61-62: Skip non-spin particles; implemented by `if ~strcmp(spin_system.comp.types{n},'S'), continue; end`.
- Lines 64-65: Field scaling factor and Zeeman frequency; implemented by `scal=(spin_system.inter.gammas(n)/gamma_chan)^2`.
- Lines 68-69: Never-resonant term for every spin; implemented by `coeff=scal/(2*(omega_n+carrier_frq(c)))`.
- Lines 71-72: Resonant-capable term for foreign isotopes; implemented by `if ~strcmp(spin_system.comp.isotopes{n},channels{c})`.
- Lines 76-77: Add the weighted longitudinal operator; implemented by `resp_op=resp_op+coeff*operator(spin_system,{'Lz'},{n})`.
- Lines 81-82: Store without clean-up, elements scale as inverse Zeeman frequencies; implemented by `resp_ops{c}=resp_op`.

### Control flow inferred from the code

- Line 50: `for` loop over `c=1:numel(channels)`.
- Line 59: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 62: conditional branch on `~strcmp(spin_system.comp.types{n},'S'), continue; end`.
- Line 72: conditional branch on `~strcmp(spin_system.comp.isotopes{n},channels{c})`.

### Key state/data transformations

- Lines 47: computes `resp_ops` using `resp_ops=cell(size(channels))`.
- Lines 53: computes `gamma_chan` using `gamma_chan=spin(channels{c})`.
- Lines 56: computes `resp_op` using `resp_op=mprealloc(spin_system,1)`.
- Lines 65: computes `scal` using `scal=(spin_system.inter.gammas(n)/gamma_chan)^2`.
- Lines 66: computes `omega_n` using `omega_n=spin_system.inter.basefrqs(n)`.
- Lines 69: computes `coeff` using `coeff=scal/(2*(omega_n+carrier_frq(c)))`.
- Lines 82: computes `resp_ops{c}` using `resp_ops{c}=resp_op`.

### Local helper functions

- Line 89: `grumble()` — `function grumble(spin_system,channels,carrier_frq)`.
  - Representative operation: `if ~isfield(spin_system,'comp')`.
  - Representative operation: `error('spin_system must be a Spinach data structure.')`.

## Parameters / inputs

- channels -cell array of isotope strings, one per
- control operator, e.g. {'1H','1H'} for
- an X,Y control operator pair
- carrier_frq -row vector of signed carrier frequencies
- (rad/s), one per control operator; for a
- transmitter on resonance this is the value
- of inter.basefrqs for the channel isotope

## Outputs

- resp_ops -cell array of Bloch-Siegert response ope-
- rators, one per control operator, in the
- formalism of the spin system provided

## Implementation structure

- Bloch-Siegert response operators for the optimal control module. For
- each control channel, returns the operator whose coefficient in every
- time slice of a GRAPE optimisation is the square of the physical con-
- trol amplitude on that channel. The operator collects the second-order
- Bloch-Siegert frequency shifts of every spin in the system:
- B=sum_n (gamma_n/gamma_c)^2*[1/(2*(omega_n+omega_c))
- +(foreign isotopes only) 1/(2*(omega_n-omega_c))]*Lz_n
- where omega_n are signed laboratory frame Zeeman frequencies from
- spin_system.inter.basefrqs and omega_c is the signed carrier frequen-
- cy of the channel. Spins belonging to the isotope that the channel
- addresses receive only the never-resonant term because their resonant
- term is the control operator itself, which GRAPE propagates exactly.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `mprealloc()`, `strcmp()`, `carrier_frq()`, `operator()`, `isfield()`, `iscell()`, `any()`, `cellfun()`, `channels()`, `all()`, `ismember()`.
