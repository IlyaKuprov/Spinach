# kernel/optimcon/bss_ops.m

- Signature: `resp_ops=bss_ops(spin_system,channels,carrier_frq)`

## Purpose

Bloch-Siegert response operators for the optimal control module. For each control channel, returns the operator whose coefficient in every time slice of a GRAPE optimisation is the square of the physical con- trol amplitude on that channel. The operator collects the second-order Bloch-Siegert frequency shifts of every spin in the system: B=sum_n (gamma_n/gamma_c)^2*[1/(2*(omega_n+omega_c)) +(foreign isotopes only) 1/

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

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
