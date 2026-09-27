# examples/fundamentals/derivative_tests/dirdiff_test_system.m

- Signature: `[spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalism)`

## Purpose

Builds the spin system and shared operators and states used by the directional-derivative tests.

## System and basis

The accepted formalisms are sphten-liouv, zeeman-liouv, and zeeman-hilb. The sphten-liouv case uses 100 non-interacting 13C spins; the two Zeeman cases use 2. The magnetic field is 28.18 T, and the chemical shifts are equally spaced over −100 to +100 ppm. For sphten-liouv, the basis uses approximation IK-2, proximity level 1, and connectivity scalar_couplings; the Zeeman formalisms use approximation none.

## Outputs

The function creates and bases the system, then returns normalized Sx, Sy, and Sz states built from the corresponding Lx, Ly, and Lz states for 13C. Each state is divided by norm(full(state),2). It also returns the 13C operators Lx and Ly, and the drift Hamiltonian hamiltonian(assume(spin_system,'nmr')).
