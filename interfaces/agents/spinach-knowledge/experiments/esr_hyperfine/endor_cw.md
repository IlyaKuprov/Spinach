# experiments/esr_hyperfine/endor_cw.m

- MATLAB implementation: [experiments/esr_hyperfine/endor_cw.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_hyperfine/endor_cw.m)

Signature: fid=endor_cw(spin_system,parameters,H,R,K)

## Purpose and physical meaning

This routine is a fast approximate simulation of isotropic continuous-wave ENDOR. It models nuclear-spin response weighted by electron–nuclear hyperfine couplings; it returns a nuclear free-induction decay (FID), whose Fourier transform approximates a CW ENDOR spectrum. It is not a magnetic-field sweep: the supplied sweep is a nuclear-frequency sweep width. The routine does not model DNP, hyperpolarisation, or spatial imaging.

## Inputs and required settings

- spin_system supplies the spin specification and is used to construct states and operators.
- parameters.sweep: nuclear frequency sweep width, in Hz. The sampled evolution step is 1/parameters.sweep seconds.
- parameters.npoints: scalar number of FID points.
- H, R, and K: Hamiltonian, relaxation, and kinetics matrices supplied by the context function. They must be matrices of identical dimensions. The function accepts the sphten-liouv and zeeman-liouv formalisms and converts to the adjoint representation when needed.

## State preparation, propagation, and signal

The routine forms L = H + 1i*R + 1i*K. For each electron–nucleus pair it computes the isotropic coupling amplitude from one third of the traces of the two coupling-matrix blocks, sums nuclear Lz states weighted by the absolute amplitudes, and normalises the resulting initial state. A nuclear Sy pulse of pi/2 radians excites the nuclei. Detection uses the nuclear L+ state. Propagation uses a dwell time of 1/sweep and npoints-1 further points; the routine returns the real part of the FID for frequency symmetrisation.

The numerical constants in the implementation are the pi/2 nuclear rotation and the sampling relation dt=1/sweep; the source does not provide a worked numeric sweep-width example. The returned value is simulated, not measured.

## Scope and limitations

This is an isotropic, approximate CW-ENDOR model, not a full field-dependent EPR acquisition. It returns an FID rather than performing the Fourier transform itself. No DOI or experimental measurement is specified by the source or the earlier page.

Source: https://spindynamics.org/wiki/index.php?title=endor_cw.m

Both the normalised hyperfine-weighted nuclear operator sum and the receiver use `coil_state`: their geometry is independent of substance concentrations, including zero populations. This approximate, normalised spectrum is not an absolute concentration-weighted signal.
