# examples/shaped_pulses/shaped_pulse_fp.m

- Signature: `shaped_pulse_fp()`

## Purpose

This example uses the Fokker–Planck shaped-pulse propagator to simulate an off-resonance rectangular soft pulse. Its source comment notes that the frequency offset accumulates phase during the pulse, as it does during a chirp. The subsequent spectrum is a simulated acquisition result, not an experimental measurement.

## Spin system and pulse

The model is a 31-proton chain at 14.1 T, with scalar Zeeman shifts from −4 to +4 ppm and 10 Hz scalar couplings between adjacent spins. It uses the IK-2 basis with scalar-coupling connectivity and proximity level 1. The initial state is 1H Lz; the controls are Lx and Ly, and the background Hamiltonian is constructed under the NMR assumption.

The call `shaped_pulse_af(spin_system,H,Lx,Ly,rho,1922.4,50.0,5e-3,-pi/2,2)` supplies a 1922.4 Hz RF frequency offset, 50.0 rad/s RF amplitude, a 5 ms pulse duration, an initial RF phase of −π/2 rad, and maximum Fokker–Planck rank 2. These scalar pulse parameters describe a single constant-amplitude, constant-frequency slice; this is not a sampled amplitude/phase waveform. The helper documentation identifies the Fokker–Planck construction with Eq. 33 of the paper at DOI [10.1016/j.jmr.2016.07.005](https://doi.org/10.1016/j.jmr.2016.07.005).

## Acquisition and observable

The pulse-prepared state is passed to liquid-state NMR acquisition with a 7000 Hz sweep, 2048 acquired points, zero filling to 8192 points, and a Hz axis. The FID is phase-corrected by exp(−i × 0.67), exponentially apodised with parameter 6, Fourier transformed, and plotted as a real spectrum. The script does not explicitly construct a relaxation superoperator or apply a gradient or homospoil step before this acquisition.

Sources: [shaped_pulse_fp.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/shaped_pulses/shaped_pulse_fp.m) and [shaped_pulse_af.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/shaped_pulse_af.m).
