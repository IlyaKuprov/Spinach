# examples/shaped_pulses/shaped_pulse_chirp_af.m

- Signature: `shaped_pulse_chirp_af()`

## Purpose

Simulate a chirp pulse in amplitude-frequency coordinates using the Fokker-Planck formalism. This method requires fewer points than the Nyquist-Shannon condition of two points per period of the highest frequency for a `{Cx,Cy}`-parameterised chirp simulation. Calculation time: seconds.

## Physical / mathematical content

- The system contains 31 `1H` spins at a 14.1 T magnetic field. Their scalar Zeeman values span `-4` to `4`, and adjacent spins have scalar couplings of `20`.
- A WURST chirp waveform is generated with `chirp_pulse(100,0.1,2000,16,'wurst')`; its frequency coordinates are then shifted by `1000`.
- The amplitude-frequency pulse acts on initial longitudinal magnetisation with `shaped_pulse_af(...,frqs,amps,durs,pi/2,2)`. A homospoil destroys stray transverse magnetisation before a hard `pi/2` pulse about `Ly` prepares detection.

## Numerical / algorithmic content

- The basis uses `sphten-liouv` formalism, `IK-2` approximation, `scalar_couplings` connectivity, and proximity level `1`.
- Acquisition detects `L+` on `1H` with zero offset, sweep `5100`, `2048` points, `8192`-point zero filling, and an axis in Hz.
- The acquired FID receives exponential apodisation with parameter `6`; `fftshift(fft(fid,parameters.zerofill))` produces the spectrum, whose real part is plotted.

## Implementation structure

- Create the spin system and basis, then construct the Hamiltonian, relaxation and kinetics operators, and the `Lx` and `Ly` operators.
- Generate and frequency-shift the chirp; apply the shaped pulse, homospoil, and hard pulse in sequence.
- Acquire, apodise, Fourier-transform, and plot the signal.
