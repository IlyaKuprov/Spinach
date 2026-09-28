# experiments/zulf/zulf_abrupt.m

- Signature: `fid=zulf_abrupt(spin_system,parameters,H,R,K)`

## Purpose

Simulates zero-field magnetometry by propagating the equilibrium initial condition through an exponential drop of the external magnetic field, then acquiring a free induction decay at zero field.

## Physical / mathematical content

The function constructs its own Zeeman and coupling Hamiltonians in the lab frame. It normalises the Zeeman Hamiltonian to 1 T, then propagates the equilibrium state through the field-drop profile using `drop(n)*H_z + H_c + 1i*R + 1i*K`. At zero field, acquisition uses `H_c + 1i*R + 1i*K`. The pulse and detection states are weighted by the spin magnetogyric ratios relative to 1H.

## Numerical / algorithmic content

The exponential field-drop profile is generated from the system magnet field, `parameters.drop_field`, `parameters.drop_time`, `parameters.drop_npoints`, and `parameters.drop_rate`; each profile value is propagated for `drop_time/drop_npoints`. The weighted transverse pulse is applied after the drop. Acquisition uses a time step of `1/parameters.sweep` and `parameters.npoints-1` observable-mode propagation steps. The detection mode is either `'quadrature'`, which retains the complex signal, or `'uniaxial'`, which returns its real part. Consistency checks require positive sweep width, acquisition point count, drop time, drop rate, and drop-step count; `drop_field` must be non-negative and flip angle a real scalar.

## Parameters / inputs

- `.drop_field` - the magnetic field that the sample should be dropped to, starting from the field specified in sys.magnet, Tesla
- `.drop_time` - drop time, seconds
- `.drop_npoints` - number of discretisation points in the drop
- `.drop_rate` - field drop rate, Hz
- `.sweep` - sweep width during acquisition
- `.npoints` - number of points during acquisition
- `.detection` - 'uniaxial' to emulate common ZULF hardware, 'quadrature' for proper frequency sign discrimination
- `.flip_angle` - pulse flip angle in radians for protons; for other nuclei, this will be scaled by the gamma ratio
- `H` - checked as a matrix argument but not used to construct the Hamiltonian; this function makes its own Hamiltonian
- `R` - relaxation superoperator, used during field drop and acquisition
- `K` - kinetics superoperator, used during field drop and acquisition

## Outputs

- `fid` - the free induction decay detected on the internally generated gamma-weighted state

Note: this function ignores the offset parameter and makes its own Hamiltonian.

## Implementation structure

After validation, the routine builds the internal Zeeman and coupling Hamiltonians, propagates the equilibrium state through the exponential field drop, constructs the gamma-weighted detection coil and pulse operator, applies the pulse, and evolves the signal at zero field. The final detection-mode switch either leaves the quadrature signal unchanged or takes its real part for uniaxial detection.
