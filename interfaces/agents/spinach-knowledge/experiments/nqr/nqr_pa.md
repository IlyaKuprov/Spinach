# experiments/nqr/nqr_pa.m

MATLAB source: [experiments/nqr/nqr_pa.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nqr/nqr_pa.m)

This is a nuclear-quadrupole-resonance soft-pulse/acquire sequence. Acquisition is described as idealised with infinite bandwidth; the pulse is a shaped, off-resonance pulse, not a DANTE train.

## Inputs and source-defined sequence

`spectrum=nqr_pa(spin_system,parameters,H,R,K)` requires `sweep` (two frequency-window limits in Hz, in ascending order), scalar positive-integer `npoints`, `rho0`, `coil`, RF operators `Lx` and `Ly`, `rf_frq` (Hz), `rf_pwr` (the multiplier in rad/s for the RF Hamiltonian), `rf_dur` (seconds), and `spc_dim`. `H`, `R`, and `K` are context-supplied matrices of matching dimensions. For a `zeeman-hilb` input, the function converts the system to Liouville space and converts `Lx` and `Ly` to commutation superoperators; it then extends them using `spc_dim`.

The pulse call uses `H+1i*R+1i*K`, both RF operators, `rho0`, and the specified frequency, amplitude, and duration. It passes the resulting state to `slowpass` for the requested frequency-domain acquisition. The output is a spectrum with `npoints` samples in the supplied window; the function does not itself construct a field-orientation or MAS sweep.

## Relaxation requirement

The source warns that relaxation must be present in the dynamics for the matrix inversion in acquisition to converge, and that `R` should not be thermalised. These are source-level input requirements, not a statement that a particular model or result has been validated.

https://spindynamics.org/wiki/index.php?title=nqr_pa.m
