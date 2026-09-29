# experiments/imaging/phase_enc_3d.m

Signature: `fid=phase_enc_3d(spin_system,parameters,H,R,K,G,F)`

## Contract and acquisition

This is an imaging callback, normally invoked as imaging(spin_system,@phase_enc_3d,parameters); the imaging framework supplies the Hamiltonian, relaxation, diffusion, gradient operators, and initial state. It is available in the sphten-liouv and zeeman-liouv formalisms. H, R, K, and F must be same-sized matrices, and G must contain at least three gradient operators.

## Sequence and spatial encoding

The background generator is B = H + F + iR + iK. G{1} is the slice-select axis: a shaped RF pulse acts through the 1H Lx and Ly operators while the slice gradient is on, followed by evolution with the opposite gradient to roll it back. The sequence evolves for t_echo, applies a 180-degree y rotation, then evolves for a second t_echo. Projection onto 1H L+ produces an internal slice profile for plotting; that intermediate is not the returned image.

The returned signal is acquired in a two-dimensional phase-encode/readout loop. For each of image_size(1) linearly spaced amplitudes from -pe_grad_amp to +pe_grad_amp, the sequence evolves with that amplitude on G{2} for pe_grad_dur. G{3} supplies the readout gradient: it first runs at negative amplitude for half ro_grad_dur, then the coil is observed under the positive readout gradient. Thus the input sample grid is three-dimensional for slice selection, while the acquired k-space array is two-dimensional.

## Required parameters, units, and output

Required fields include rho0 (initial state), coil (numeric detection state), npts (positive-integer spatial grid), dims (three real sample extents), image_size (two odd integers, each at least three; `imaging()` rejects even sizes), ss_grad_amp, pe_grad_amp, ro_grad_amp (gradient amplitudes in T/m), pe_grad_dur, ro_grad_dur, t_echo (times in seconds), rf_frq_list (Hz), rf_amp_list (rad/s), rf_dur_list (seconds), and rf_phi (pulse phase; the example uses pi/2). The RF frequency, amplitude, and duration vectors must have matching lengths. image_size(1) is the number of phase-encode lines; image_size(2) is the number of readout samples. The readout dwell is ro_grad_dur/(image_size(2)-1) seconds. The function returns complex fid with shape image_size, containing raw k-space samples, not a reconstructed or measured image.

## Source-backed configuration example

The supplied high-resolution example uses image_size [129 129], 32 mT/m for slice, phase-encode, and readout gradients, pe_grad_dur 0.2 ms, ro_grad_dur 0.3 ms, and t_echo 20 ms. Its slice RF table has 50 intervals totaling 0.2 ms and uses a Gaussian amplitude shape. These are example inputs in examples/extremes/ph_enc_3d_highres.m; no simulation or measured image is claimed here.

## References

- [Spinach Wiki: phase_enc_3d.m](https://spindynamics.org/wiki/index.php?title=phase_enc_3d.m)
- [Canonical MATLAB source: experiments/imaging/phase_enc_3d.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/phase_enc_3d.m)
- [Example configuration: ph_enc_3d_highres.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/extremes/ph_enc_3d_highres.m)
