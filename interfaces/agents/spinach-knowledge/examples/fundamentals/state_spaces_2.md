# examples/fundamentals/state_spaces_2.m

- Signature: `state_spaces_2()`

## Purpose

Tracks spin-correlation-order contributions during a pulse-acquire `1H` NMR simulation of anti-3,5-difluoroheptane. The source describes a 16-spin example and reports that correlations through order eight suffice for practical simulation in this case; nine- and ten-spin contributions are also shown.

## Physical / mathematical content

- The model uses the source's 11.7464 T spin-system parameters, including the listed `12C`, `1H`, and `19F` isotopes, chemical shifts, and scalar couplings.
- The manually restricted IK-0 basis uses inter-level 1, specified spin subsets and symmetries, longitudinal `19F`, and projection +1. Automatic state dropout is disabled with `sys.disable={'zte'}`.

## Numerical / algorithmic content

- Starts from `L+` on `1H`, propagates the trajectory for 1000 steps at 1 ms with the NMR-assumed Hamiltonian, and analyses it with `trajan(...,'correlation_order')`.
- The source estimates minutes of runtime and notes that a GPU is faster; its GPU-enable line is commented out.

## Implementation structure

- Defines the spin system and manually specified basis, computes the trajectory, and plots the correlation-order contributions on a logarithmic scale.
