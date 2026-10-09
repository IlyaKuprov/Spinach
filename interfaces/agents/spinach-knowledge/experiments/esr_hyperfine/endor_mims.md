# experiments/esr_hyperfine/endor_mims.m

- MATLAB implementation: [experiments/esr_hyperfine/endor_mims.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_hyperfine/endor_mims.m)

Signature: fid=endor_mims(spin_system,parameters,H,R,K)

## Purpose and physical sequence

This routine simulates Mims ENDOR using ideal hard pulses. It prepares electron Lz magnetisation, applies electron pi/2 – tau – pi/2 pulses, selects zero electron coherence, applies a nuclear pi/2 pulse, and selects nuclear coherence orders -1 and +1. It propagates the indirect nuclear-frequency dimension, applies the phase-difference nuclear-pulse operation, then applies an electron pi/2 refocusing pulse and detects the stimulated echo on the electron L+ state after tau. The returned FID Fourier-transforms to a Mims ENDOR signal.

This is a simulated hyperfine-sensitive electron–nuclear sequence, not a measured acquisition. Its sweep parameter is a nuclear-frequency sweep width; the function does not sweep the static magnetic field and does not implement DNP, hyperpolarisation, or spatial imaging.

## Inputs and numerical settings

- spin_system; H, R, and K are the spin system and context-provided Hamiltonian, relaxation, and kinetics matrices. The matrices must have identical dimensions. Supported formalisms are sphten-liouv and zeeman-liouv; conversion to the adjoint representation is performed if needed.
- parameters.sweep: nuclear frequency sweep width, in Hz. It sets the indirect-dimension time step to 1/sweep seconds.
- parameters.npoints: number of FID points.
- parameters.tau: stimulated-echo time, in seconds.

The implementation uses pi/2 rotations and propagates npoints-1 further indirect-dimension points. The source gives no worked numeric sweep-width or tau example. Relaxation and kinetics enter through L = H + 1i*R + 1i*K; the output is the electron-coil-detected FID, not a spectrum already Fourier transformed.

Source: https://spindynamics.org/wiki/index.php?title=endor_mims.m
