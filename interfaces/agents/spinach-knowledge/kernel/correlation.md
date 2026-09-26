# kernel/correlation.m

- Signature: `rho=correlation(spin_system,rho,orders,spins)`

## Purpose

Correlation order selection function -keeps only the specified orders of spin correlation in the state vector. This is useful as an analyti- cal replacement for complicated phase cycles. Syntax: rho=correlation(spin_system,rho,correlation_orders,spins)

## Physical / mathematical content

- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

## Parameters / inputs

- rho -a state vector or a horizontal stack thereof;
- in zeeman-hilb, a density matrix or a horizon-
- tal stack thereof
- correlation_orders -a row vector of correlation
- orders to keep
- spins -which spins to consider (e.g.
- '1H', '13C', 'all')

## Outputs

- rho -the state vector with the undesired orders of
- spin correlations zeroed out
- Note: this function requires sphten-liouv, zeeman-liouv, or zeeman-
- hilb formalism; Fokker-Planck direct products are supported in
- the Liouville space formalisms. In the Zeeman formalisms the
- selection is an exact projection built from per-spin identity
- component channels because correlation order is not diagonal
- in the Zeeman basis; zeeman-hilb density matrices are stretch-
- ed into Liouville space, filtered there, and folded back.

## Implementation structure

- Correlation order selection function -keeps only the specified orders
- of spin correlation in the state vector. This is useful as an analyti-
- cal replacement for complicated phase cycles. Syntax:
- rho=correlation(spin_system,rho,correlation_orders,spins)
- rho - a state vector or a horizontal stack thereof;
- in zeeman-hilb, a density matrix or a horizon-
- tal stack thereof
- correlation_orders - a row vector of correlation
- orders to keep
- spins - which spins to consider (e.g.
- '1H', '13C', 'all')
- rho -the state vector with the undesired orders of
