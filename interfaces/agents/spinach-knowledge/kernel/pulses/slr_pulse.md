# kernel/pulses/slr_pulse.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/slr_pulse.m`
- Signature: `[Cx,Cy,durs,amps,phis]=slr_pulse(npts,dur,tbw,flip_angle,pass_rip,stop_rip)`
- Total lines: 193

## Purpose

Shinnar-Le Roux linear-phase selective excitation pulse. Syntax: [Cx,Cy,durs,amps,phis]=slr_pulse(npts,dur,tbw,flip_angle,pass_rip,stop_rip)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- npts -even number of piecewise-constant pulse slices
- dur -total pulse duration, seconds
- tbw -time-bandwidth product, defined as the pulse
- duration times the nominal full passband width
- flip_angle -on-resonance flip angle between zero and pi/2,
- radians
- pass_rip -90-degree excitation passband ripple target used
- in the prototype design, dimensionless
- stop_rip -90-degree excitation stopband ripple target used
- in the prototype design, dimensionless

## Outputs

- Cx -X control amplitudes, rad/s, 1 x npts row vector
- Cy -Y control amplitudes, rad/s, 1 x npts row vector
- durs -pulse slice durations, seconds, 1 x npts row vector
- amps -RF amplitudes, rad/s, 1 x npts row vector
- phis -RF phases, radians, 1 x npts row vector
- The beta polynomial is obtained by continuous weighted least squares
- in a linear-phase cosine basis. The complementary minimum-phase alpha
- polynomial and the RF waveform are then obtained by the inverse SLR
- transform. The ripple arguments enter the excitation-pulse transform
- and the transition-width estimate of Pauly et al.; they are design
- targets rather than guaranteed minimax error bounds.
- For flip angles below pi/2, they do not specify angle-independent
- magnetisation error bounds.
- The output controls are calibrated for Spinach propagation under
- exp(-1i*H*t) and may be passed directly to shaped_pulse_xy().
- J. Pauly, P. Le Roux, D. Nishimura, and A. Macovski,
- IEEE Transactions on Medical Imaging 10(1), 53-65 (1991),

## Implementation structure

- Shinnar-Le Roux linear-phase selective excitation pulse. Syntax:
- [Cx,Cy,durs,amps,phis]=slr_pulse(npts,dur,tbw,flip_angle,pass_rip,stop_rip)
- npts -even number of piecewise-constant pulse slices
- dur -total pulse duration, seconds
- tbw -time-bandwidth product, defined as the pulse
- duration times the nominal full passband width
- flip_angle -on-resonance flip angle between zero and pi/2,
- radians
- pass_rip -90-degree excitation passband ripple target used
- in the prototype design, dimensionless
- stop_rip -90-degree excitation stopband ripple target used
- Cx -X control amplitudes, rad/s, 1 x npts row vector

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `log10()`, `integral()`, `basis_fun()`, `decomposition()`, `flipud()`, `eps()`, `nextpow2()`, `beta_pad()`, `any()`, `log_spec()`, `fliplr()`, `alpha_poly()`, `alpha()`, `beta()`, `hypot()`.
