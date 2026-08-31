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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 57-58: Check consistency; implemented by `grumble(npts,dur,tbw,flip_angle,pass_rip,stop_rip)`.
- Lines 60-61: Convert 90-degree excitation-profile targets into beta targets; implemented by `beta_pass=sqrt(pass_rip/2)`.
- Lines 65-66: Estimate the transition width using the Pauly relation; implemented by `log_pass=log10(beta_pass)`.
- Lines 73-74: Reject filter specifications that do not fit below Nyquist; implemented by `if (pass_edge<=0)||(pass_edge>=stop_edge)||(stop_edge>=1)`.
- Lines 78-79: Build the continuous least-squares system; implemented by `orders=(0:(npts/2-1))+0.5`.
- Lines 90-91: Solve for the linear-phase beta polynomial; implemented by `cos_coeff=decomposition((gram+gram')/2,'chol')\moment`.
- Lines 95-96: Set the beta-polynomial centre response to the requested flip; implemented by `centre_resp=sum(beta)`.
- Lines 102-103: Sample the beta response on an oversampled Fourier grid; implemented by `fft_len=2^(nextpow2(npts)+4)`.
- Lines 109-110: Require a strictly positive complementary spectrum; implemented by `if resp_peak>=1`.
- Lines 118-119: Obtain the minimum-phase complementary response by cepstral factorisation; implemented by `alpha_mag=sqrt(complement)`.
- Lines 127-128: Remove one Cayley-Klein rotation at a time; implemented by `rf_angles=zeros(1,npts)`.
- Lines 150-151: Convert per-slice rotations into Spinach control amplitudes; implemented by `slice_dur=dur/npts`.
- Lines 157-158: Return Cartesian, timing, and polar waveform coordinates; implemented by `Cx=real(rf_ctrl)`.

### Control flow inferred from the code

- Line 74: conditional branch on `(pass_edge<=0)||(pass_edge>=stop_edge)||(stop_edge>=1)`.
- Line 97: conditional branch on `abs(centre_resp)<=eps(norm(beta,1))`.
- Line 110: conditional branch on `resp_peak>=1`.
- Line 114: conditional branch on `any(complement<=0)`.
- Line 129: `for` loop over `slice=npts:-1:1`.
- Line 130: conditional branch on `abs(alpha(slice))<=eps(norm(alpha,inf))`.
- Line 137: conditional branch on `sin_size==0`.
- Line 142: conditional branch on `slice>1`.
- Line 153: conditional branch on `any(~isfinite(rf_ctrl))`.

### Key state/data transformations

- Lines 61: computes `beta_pass` using `beta_pass=sqrt(pass_rip/2)`.
- Lines 62: computes `beta_stop` using `beta_stop=stop_rip/sqrt(2)`.
- Lines 63: computes `beta_scale` using `beta_scale=sin(flip_angle/2)`.
- Lines 66: computes `log_pass` using `log_pass=log10(beta_pass)`.
- Lines 67: computes `log_stop` using `log_stop=log10(beta_stop)`.
- Lines 68-69: computes `trans_measure` using `trans_measure=(5.309e-3*log_pass^2+7.114e-2*log_pass-4.761e-1)*log_stop+ (-2.66e-3*log_pass^2-5.941e-1*log_pass-4.278e-1)`.
- Lines 70: computes `pass_edge` using `pass_edge=(tbw-trans_measure)/npts`.
- Lines 71: computes `stop_edge` using `stop_edge=(tbw+trans_measure)/npts`.
- Lines 79: computes `orders` using `orders=(0:(npts/2-1))+0.5`.
- Lines 80: computes `pass_lim` using `pass_lim=pi*pass_edge`.
- Lines 81: computes `stop_lim` using `stop_lim=pi*stop_edge`.
- Lines 82: computes `basis_fun` using `basis_fun=@(x)cos(x*orders)`.
- Lines 83: computes `stop_weight` using `stop_weight=beta_pass/beta_stop`.
- Lines 84-87: computes `gram` using `gram=integral(@(x)basis_fun(x)'*basis_fun(x),0,pass_lim, 'ArrayValued',true)+stop_weight* integral(@(x)basis_fun(x)'*basis_fun(x),stop_lim,pi, 'ArrayValued',true)`.
- Lines 88: computes `moment` using `moment=integral(@(x)basis_fun(x)',0,pass_lim,'ArrayValued',true)`.
- Lines 91: computes `cos_coeff` using `cos_coeff=decomposition((gram+gram')/2,'chol')\moment`.
- Lines 92: computes `left_half` using `left_half=flipud(cos_coeff)/2`.
- Lines 93: computes `beta` using `beta=[left_half;flipud(left_half)]'`.

### Local helper functions

- Line 167: `grumble()` — `function grumble(npts,dur,tbw,flip_angle,pass_rip,stop_rip)`.
  - Representative operation: `if (~isnumeric(npts))||(~isreal(npts))||(~isscalar(npts))|| (~isfinite(npts))||(npts<2)||(mod(npts,2)~=0)`.
  - Representative operation: `(~isfinite(npts))||(npts<2)||(mod(npts,2)~=0)`.

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
