# experiments/echo_sweep.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/echo_sweep.m`
- Signature: `echo=echo_sweep(spin_system,parameters,H,~,~)`
- Total lines: 301

## Purpose

Two-pulse echo-detected frequency-swept experiment, static or under magic angle spinning, in Hilbert space, written for the EPR case of a spinning P1 centre in diamond. Two pulses of equal duration are separated by a delay, the carrier is stepped across the sweep, and the complex echo integral is returned at each carrier offset. The sequence steps through the Hamiltonian rotor stack supplied by `singlerot()`: at each time step, the stack element nearest to the rotor phase at the middle of the step is used, and the rotor phase at the start of the sequence, which stands in for the crystallite azimuth about the rotor axis, is averaged over. The coherence pathway of the pulsed spin (-1 after the first pulse, +1 after the second) is selected in place of a phase cycle.

## Physical / mathematical content

- Pulsed EPR under sample spinning: anisotropic hyperfine and Zeeman interactions make the resonance frequency of each spin packet move during the sequence, so lines with anisotropy dephase between the pulses and the echo whereas isotropic lines survive.
- The carrier offset enters as `2*pi*offset*Lz` on the pulsed spin; the free-evolution propagators are computed once per orientation and the offset is applied as a separate propagator, which is exact only when the rotor stack commutes with the `Lz` of that spin, as it does for an electron under the `'esr'` assumption set of the context function.
- The rotor phase at the start of the sequence is equivalent to the crystallite azimuth about the rotor axis, so averaging over `parameters.nphases` start positions spread evenly through the stack lets a two-angle powder grid be used.

## Numerical / algorithmic content

- Time propagation is explicit in Hilbert space: one `propagator()` per rotor stack element for free evolution, one per element and per carrier offset for the pulses, each applied as `P*rho*P'` at every time step of `parameters.timestep`.
- The stack index at each step is `round(rate*timestep*(s-1/2)*spc_dim)` past the start index, taken modulo `spc_dim`, so `parameters.rate=0` is the static case with no special code; the stack is a table of the Hamiltonian against rotor phase whose resolution should match the time step at the fastest rate used, `parameters.max_rank` of the context of about `1/(2*abs(rate)*timestep)` at that rate; a finer stack costs propagators for no gain beyond the time step, a coarser one loses phase resolution, and at slower rates consecutive steps reuse a stack element.
- The coherence pathway is enforced with `coherence()` after each pulse; the echo is the sum of `trace(coil'*rho)` over the `echo_win` steps after the second pulse, multiplied by `parameters.timestep` (Riemann sum) and divided by `parameters.nphases` at the end.
- The file contains an explicit `grumble(...)` validator that checks the formalism, the stack, the states, every sequence parameter, and the commutation of the stack with the `Lz` of the pulsed spin against `spin_system.tols.liouv_zero`.

## Parameters / inputs

- `parameters.spins` - one-element cell array naming the spin the pulses are applied to, e.g. `{'E'}`
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
- `H` - vector cell array of Hamiltonian matrices, one for each rotor phase, received from context function
- `R`, `K` - relaxation and kinetics superoperators, received from context function, not used

## Outputs

- `echo` - complex echo signal integrated over the echo window (sum over the time steps multiplied by the time step) and averaged over the rotor phases at the start of the sequence, at each carrier offset, a column vector with `parameters.npoints` elements

## Notes

- The elements of the rotor stack must commute with the `Lz` operator of the pulsed spin, as they do for an electron under the `'esr'` assumption set, because the carrier offset is applied as a separate propagator and only the pulse propagators are rebuilt at each carrier offset.
- The rotor stack is a table of the Hamiltonian against the rotor phase; its resolution should match the time step at the fastest spinning rate used, `parameters.max_rank` of the context function of about `1/(2*abs(rate)*timestep)` at that rate. A finer stack costs propagators without gaining accuracy beyond the time step, a coarser one loses rotor phase resolution; at slower rates consecutive steps reuse elements.
- The sequence is not restricted to electrons: `parameters.spins` names the spin that is pulsed and whose coherence pathway is selected, so a nuclear two-pulse echo is obtained by naming the nucleus.
- Used by `examples/esr_sol_pulsed/mas_diamond_p1.m` through `singlerot()` in the `zeeman-hilb` formalism.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ft_axis()`, `operator()`, `propagator()`, `coherence()`, `trace()`.
