# experiments/echo_sweep.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/echo_sweep.m`
- Signature: `echo=echo_sweep(spin_system,parameters,H,~,~)`
- Total lines: 288

## Purpose

Two-pulse echo-detected frequency-swept EPR experiment, static or under magic angle spinning, in Hilbert space. Two pulses of equal duration are separated by a delay, the carrier is stepped across the sweep, and the complex echo integral is returned at each carrier offset. The sequence steps through the Hamiltonian rotor stack supplied by `singlerot()`: at each time step, the stack element nearest to the rotor phase at the middle of the step is used, and the rotor phase at the start of the sequence, which stands in for the crystallite azimuth about the rotor axis, is averaged over. The electron coherence pathway (-1 after the first pulse, +1 after the second) is selected in place of a phase cycle.

## Physical / mathematical content

- Pulsed EPR under sample spinning: anisotropic hyperfine and Zeeman interactions make the resonance frequency of each spin packet move during the sequence, so lines with anisotropy dephase between the pulses and the echo whereas isotropic lines survive.
- The carrier offset enters as `2*pi*offset*Lz` on the electron; the free-evolution propagators are computed once per orientation and the offset is applied as a separate propagator, which is exact only when the rotor stack commutes with the electron `Lz`, as it does under the `'esr'` assumption set of the context function.
- The rotor phase at the start of the sequence is equivalent to the crystallite azimuth about the rotor axis, so averaging over `parameters.nphases` start positions spread evenly through the stack lets a two-angle powder grid be used.

## Numerical / algorithmic content

- Time propagation is explicit in Hilbert space: one `propagator()` per rotor stack element for free evolution, one per element and per carrier offset for the pulses, each applied as `P*rho*P'` at every time step of `parameters.timestep`.
- The stack index at each step is `round(rate*timestep*(s-1/2)*spc_dim)` past the start index, taken modulo `spc_dim`, so `parameters.rate=0` is the static case with no special code; the stack should advance by at most one element per time step at the fastest rate used, so that no rotor phase is skipped, which asks for `parameters.max_rank` of the context of at least `1/(2*rate*timestep)` at that rate (at slower rates consecutive steps reuse a stack element).
- The coherence pathway is enforced with `coherence()` after each pulse; the echo is the sum of `trace(coil'*rho)` over the `echo_win` steps after the second pulse, accumulated over the start phases and divided by `parameters.nphases` at the end.
- The file contains an explicit `grumble(...)` validator that checks the formalism, the stack, the states, every sequence parameter, and the commutation of the stack with the electron `Lz` against `spin_system.tols.liouv_zero`.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 86-87: Check consistency; implemented by `grumble(spin_system,parameters,H)`.
- Lines 89-90: Carrier offsets across the sweep; implemented by `offsets=ft_axis(0,parameters.sweep,parameters.npoints)`.
- Lines 92-93: Pulse and offset operators; implemented by `sx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 96-97: Step counts of the pulses, the delay, and the echo window; implemented by `pulse_steps=round(parameters.pulse_dur/parameters.timestep)`.
- Lines 102-103: Rotor stack advance at the middle of each time step; implemented by `stack_shift=round(parameters.rate*parameters.timestep*((1:nsteps)-1/2)*parameters.spc_dim)`.
- Lines 105-106: Rotor stack indices at the start of the sequence; implemented by `start_idx=floor((0:(parameters.nphases-1))*parameters.spc_dim/parameters.nphases)`.
- Lines 108-109: Free evolution propagators at every rotor phase; implemented by `p_free=cell(parameters.spc_dim,1)`.
- Lines 114-115: Preallocate the answer; implemented by `echo=zeros(parameters.npoints,1)`.
- Lines 117-118: Loop over carrier offsets; implemented by `for k=1:parameters.npoints`.
- Lines 120-121: Carrier offset propagator; implemented by `p_off=propagator(spin_system,2*pi*offsets(k)*sz,parameters.timestep)`.
- Lines 123-124: Pulse propagators at every rotor phase; implemented by `p_pulse=cell(parameters.spc_dim,1)`.
- Lines 130-131: Loop over rotor phases at the start of the sequence; implemented by `for j=1:parameters.nphases`.
- Lines 133-134: Rotor stack indices at each time step; implemented by `idx=mod(start_idx(j)+stack_shift,parameters.spc_dim)+1`.
- Lines 136-137: First pulse; implemented by `rho=parameters.rho0`.
- Lines 142-143: Select the -1 coherence on the electron; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},-1}})`.
- Lines 145-146: Interpulse delay; implemented by `for s=(pulse_steps+1):(pulse_steps+delay_steps)`.
- Lines 150-151: Second pulse; implemented by `for s=(pulse_steps+delay_steps+1):(2*pulse_steps+delay_steps)`.
- Lines 155-156: Select the +1 coherence on the electron; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},+1}})`.
- Lines 158-159: Integrate the signal over the echo window; implemented by `for s=(2*pulse_steps+delay_steps+1):nsteps`.
- Lines 168-169: Average over the rotor phases at the start of the sequence; implemented by `echo=echo/parameters.nphases`.

### Control flow inferred from the code

- Line 110: `for` loop over `n=1:parameters.spc_dim`.
- Line 118: `for` loop over `k=1:parameters.npoints`.
- Line 125: `for` loop over `n=1:parameters.spc_dim`.
- Line 131: `for` loop over `j=1:parameters.nphases`.
- Line 138: `for` loop over `s=1:pulse_steps`.
- Line 146: `for` loop over `s=(pulse_steps+1):(pulse_steps+delay_steps)`.
- Line 151: `for` loop over `s=(pulse_steps+delay_steps+1):(2*pulse_steps+delay_steps)`.
- Line 159: `for` loop over `s=(2*pulse_steps+delay_steps+1):nsteps`.

### Key state/data transformations

- Lines 90: computes `offsets` using `offsets=ft_axis(0,parameters.sweep,parameters.npoints)`.
- Lines 93-94: computes `sx` and `sz` using `operator(spin_system,'Lx',parameters.spins{1})` and `operator(spin_system,'Lz',parameters.spins{1})`.
- Lines 97-100: computes `pulse_steps`, `delay_steps`, `echo_steps`, and `nsteps` using `round(parameters.pulse_dur/parameters.timestep)`, `round(parameters.tau/parameters.timestep)`, `round(parameters.echo_win/parameters.timestep)`, and `nsteps=2*pulse_steps+delay_steps+echo_steps`.
- Lines 103: computes `stack_shift` using `stack_shift=round(parameters.rate*parameters.timestep*((1:nsteps)-1/2)*parameters.spc_dim)`.
- Lines 106: computes `start_idx` using `start_idx=floor((0:(parameters.nphases-1))*parameters.spc_dim/parameters.nphases)`.
- Lines 111: computes `p_free{n}` using `p_free{n}=propagator(spin_system,H{n},parameters.timestep)`.
- Lines 115: computes `echo` using `echo=zeros(parameters.npoints,1)`.
- Lines 121: computes `p_off` using `p_off=propagator(spin_system,2*pi*offsets(k)*sz,parameters.timestep)`.
- Lines 126-127: computes `p_pulse{n}` using `p_pulse{n}=propagator(spin_system,H{n}+2*pi*offsets(k)*sz+ 2*pi*parameters.pulse_frq*sx,parameters.timestep)`.
- Lines 134: computes `idx` using `idx=mod(start_idx(j)+stack_shift,parameters.spc_dim)+1`.
- Lines 139: computes `rho` using `rho=p_pulse{idx(s)}*rho*p_pulse{idx(s)}'`.
- Lines 147: computes `rho` using `rho=p_off*p_free{idx(s)}*rho*p_free{idx(s)}'*p_off'`.
- Lines 161: computes `echo(k)` using `echo(k)=echo(k)+trace(parameters.coil'*rho)`.
- Lines 169: computes `echo` using `echo=echo/parameters.nphases`.

### Local helper functions

- Line 174: `grumble()` — `function grumble(spin_system,parameters,H)`. Requires the `zeeman-hilb` formalism, a cell array `H` of `parameters.spc_dim` square matrices of one dimension that commute with the electron `Lz` to within `spin_system.tols.liouv_zero`, `rho0` and `coil` of that dimension, positive real scalar `pulse_dur`, `pulse_frq`, `echo_win`, `timestep`, and `sweep` with `pulse_dur` and `echo_win` not shorter than `timestep`, non-negative real scalar `tau`, real scalar `rate`, positive integer `nphases` not exceeding `spc_dim`, and integer `npoints` greater than two as `ft_axis()` requires.

## Parameters / inputs

- `parameters.spins` - one-element cell array with the electron specification, e.g. `{'E'}`
- `parameters.rho0` - initial state, a density matrix
- `parameters.coil` - detection state, a density matrix
- `parameters.pulse_dur` - duration of each pulse, seconds
- `parameters.pulse_frq` - nutation frequency of the pulses, Hz
- `parameters.tau` - delay between the end of the first pulse and the start of the second pulse, seconds
- `parameters.echo_win` - echo integration window after the end of the second pulse, seconds
- `parameters.timestep` - propagation time step, seconds; the pulses, the delay, and the echo window are rounded to whole steps
- `parameters.rate` - spinning rate, Hz, zero for a static sample
- `parameters.nphases` - number of rotor phases at the start of the sequence to average over
- `parameters.sweep` - width of the carrier sweep, Hz
- `parameters.npoints` - number of carrier offsets, placed on the `ft_axis` grid of the sweep
- `parameters.spc_dim` - number of elements in the rotor stack, received from context function
- `H` - cell array of Hamiltonian matrices, one for each rotor phase, received from context function
- `R`, `K` - relaxation and kinetics superoperators, received from context function, not used

## Outputs

- `echo` - complex echo integral at each carrier offset, averaged over the rotor phases at the start of the sequence, a column vector with `parameters.npoints` elements

## Notes

- The elements of the rotor stack must commute with the electron `Lz` operator, as they do under the `'esr'` assumption set, because the carrier offset is applied as a separate propagator and only the pulse propagators are rebuilt at each carrier offset.
- The rotor stack should be fine enough to advance by at most one element per time step at the fastest spinning rate used, so that no rotor phase is skipped; `parameters.max_rank` of the context function should be at least `1/(2*rate*timestep)` at that rate.
- Used by `examples/esr_sol_pulsed/mas_diamond_p1.m` through `singlerot()` in the `zeeman-hilb` formalism.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ft_axis()`, `operator()`, `propagator()`, `coherence()`, `trace()`.
