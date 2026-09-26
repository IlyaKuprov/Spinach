# experiments/zulf/zulf_abrupt.m

- Signature: `fid=zulf_abrupt(spin_system,parameters,H,R,K)`

## Purpose

Zero-field magnetometry experiment that propagates the initial condition through an exponential drop in the external magnetic field and then runs the detection at zero field. Syntax: fid=zulf_abrupt(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Zero-field experiment implementations. They propagate J-coupled systems in the absence of strong carrier terms and often model abrupt field switching.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Parameters / inputs

- .drop_field -the magnetic field that the sample should
- be dropped to, starting from the field spe-
- cified in sys.magnet, Tesla
- .drop_time -drop time, seconds
- .drop_npoints -number of discretisation points in the drop
- .drop_rate -field drop rate, Hz
- .sweep -sweep width during acquisition
- .npoints -number of points during acquisition
- .detection -'uniaxial' to emulate common ZULF
- hardware, 'quadrature' for proper
- frequency sign discrimination
- .flip_angle -pulse flip angle in radians for
- protons; for other nuclei, this
- will be scaled by the gamma ratio

## Outputs

- fid -the free induction decay detected on the
- internally generated gamma-weighted state
- Note: this function ignores the offset parameter and makes its
- own Hamiltonian.

## Implementation structure

- Zero-field magnetometry experiment that propagates the initial condition
- through an exponential drop in the external magnetic field and then runs
- the detection at zero field. Syntax:
- fid=zulf_abrupt(spin_system,parameters,H,R,K)
- .drop_field -the magnetic field that the sample should
- be dropped to, starting from the field spe-
- cified in sys.magnet, Tesla
- .drop_time -drop time, seconds
- .drop_npoints -number of discretisation points in the drop
- .drop_rate -field drop rate, Hz
- .sweep -sweep width during acquisition
- .npoints -number of points during acquisition
