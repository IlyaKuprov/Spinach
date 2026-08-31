# examples/optimal_control/bloch_siegert/yusuke_14n_broadening_demo.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/bloch_siegert/yusuke_14n_broadening_demo.m`
- Signature: `yusuke_14n_broadening_demo()`
- Total lines: 205

## Purpose

Reduced effective-model illustration of the trade-off discussed by Nehra, Agarwal, and Nishiyama for 14N decoupling under 1H detection. The example is intentionally qualitative rather than quantitative: it shows why low-power CW decoupling is narrowband, why increasing the 14N RF field produces a Bloch-Siegert shift on the observed 1H resonance, and why B1 inhomogeneity turns that shift into broadening. The "offset-t

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.

## Numerical / algorithmic content

- The file also defines local helper function(s): `decoupling_profile()`, `proton_line()`, `lorentzian()`, `gaussian_weights()`, `normalise_line()`, `linewidth()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Magnetic field corresponding to 800 MHz 1H; implemented by `magnet=18.8`.
- Lines 21-22: Three inequivalent 14N sites around the decoupler carrier; implemented by `site_offsets_hz=[-18e3 0 18e3]`.
- Lines 25-26: B1 distribution across the sample volume; implemented by `b1_scales=linspace(0.85,1.15,81)`.
- Lines 29-30: Frequency axes for the plots; implemented by `offset_axis_hz=linspace(-40e3,40e3,801)`.
- Lines 33-34: Effective decoupling settings; implemented by `rf_low_hz=8e3`.
- Lines 38-39: Reduced-model parameters controlling the observed 1H line; implemented by `intrinsic_fwhm_hz=80`.
- Lines 43-44: Decoupling profiles across the 14N offset range; implemented by `eta_low_cw=decoupling_profile(offset_axis_hz,rf_low_hz,1.0,2)`.
- Lines 48-51: Predicted 1H line shapes; implemented by `[line_low,stats_low]=proton_line(freq_axis_hz,site_offsets_hz, site_weights,b1_scales,b1_weights,rf_low_hz,1.0,2, intrinsic_fwhm_hz,residual_penalty_hz,bs_coeff_hz)`.
- Lines 65-66: Plot the reduced-model story; implemented by `figure`.
- Lines 96-97: Console summary; implemented by `fprintf('\n')`.

### Control flow inferred from the code

- Line 72: `for` loop over `n=1:numel(site_offsets_hz)`.

### Key state/data transformations

- Lines 18: computes `magnet` using `magnet=18.8`.
- Lines 19: computes `nu_n14_hz` using `nu_n14_hz=abs(spin('14N')*magnet/(2*pi))`.
- Lines 22: computes `site_offsets_hz` using `site_offsets_hz=[-18e3 0 18e3]`.
- Lines 23: computes `site_weights` using `site_weights=[0.30 0.40 0.30]`.
- Lines 26: computes `b1_scales` using `b1_scales=linspace(0.85,1.15,81)`.
- Lines 27: computes `b1_weights` using `b1_weights=gaussian_weights(b1_scales,1.0,0.06)`.
- Lines 30: computes `offset_axis_hz` using `offset_axis_hz=linspace(-40e3,40e3,801)`.
- Lines 31: computes `freq_axis_hz` using `freq_axis_hz=linspace(-800,800,4001)`.
- Lines 34: computes `rf_low_hz` using `rf_low_hz=8e3`.
- Lines 35: computes `rf_high_hz` using `rf_high_hz=20e3`.
- Lines 36: computes `rf_ot_hz` using `rf_ot_hz=12e3`.
- Lines 39: computes `intrinsic_fwhm_hz` using `intrinsic_fwhm_hz=80`.
- Lines 40: computes `residual_penalty_hz` using `residual_penalty_hz=220`.
- Lines 41: computes `bs_coeff_hz` using `bs_coeff_hz=8e-7`.
- Lines 44: computes `eta_low_cw` using `eta_low_cw=decoupling_profile(offset_axis_hz,rf_low_hz,1.0,2)`.
- Lines 45: computes `eta_high_cw` using `eta_high_cw=decoupling_profile(offset_axis_hz,rf_high_hz,1.0,2)`.
- Lines 46: computes `eta_low_ot` using `eta_low_ot=decoupling_profile(offset_axis_hz,rf_ot_hz,2.2,4)`.
- Lines 49-51: computes `[line_low,stats_low]` using `[line_low,stats_low]=proton_line(freq_axis_hz,site_offsets_hz, site_weights,b1_scales,b1_weights,rf_low_hz,1.0,2, intrinsic_fwhm_hz,residual_penalty_hz,bs_coeff_hz)`.

### Local helper functions

- Line 115: `decoupling_profile()` — `function eta=decoupling_profile(offset_hz,rf_hz,bandwidth_gain,profile_order)`. Simple offset-response model: wider bandwidth at larger RF field or more sophisticated modulation, flatter passband at higher order.
  - Representative operation: `bw_hz=bandwidth_gain*rf_hz`.
  - Representative operation: `eta=1./(1+(abs(offset_hz)./bw_hz).^profile_order)`.
- Line 124: `proton_line()` — `function [line_shape,stats]=proton_line(freq_axis_hz,site_offsets_hz,`. Normalise the ensemble weights
  - Representative operation: `site_weights,b1_scales,b1_weights,rf_hz,bandwidth_gain,profile_order, intrinsic_fwhm_hz,residual_penalty_hz,bs_coeff_hz)`.
  - Representative operation: `intrinsic_fwhm_hz,residual_penalty_hz,bs_coeff_hz)`.
- Line 171: `lorentzian()` — `function y=lorentzian(x,centre_hz,fwhm_hz)`. Lorentzian line with unit integral
  - Representative operation: `gamma=fwhm_hz/2`.
  - Representative operation: `y=(gamma/pi)./((x-centre_hz).^2+gamma^2)`.
- Line 179: `gaussian_weights()` — `function weights=gaussian_weights(grid,centre,width)`. Truncated Gaussian weights for the B1 distribution
  - Representative operation: `weights=exp(-((grid-centre)/width).^2/2)`.
  - Representative operation: `weights=weights/sum(weights)`.
- Line 187: `normalise_line()` — `function line_out=normalise_line(line_in)`. Normalise for plotting
  - Representative operation: `line_out=line_in/max(line_in)`.
- Line 194: `linewidth()` — `function fwhm_hz=linewidth(freq_axis_hz,line_shape)`. Numerical full width at half maximum
  - Representative operation: `line_shape=line_shape/max(line_shape)`.
  - Representative operation: `idx=find(line_shape>=0.5)`.

## Implementation structure

- Reduced effective-model illustration of the trade-off discussed by
- Nehra, Agarwal, and Nishiyama for 14N decoupling under 1H detection.
- The example is intentionally qualitative rather than quantitative:
- it shows why low-power CW decoupling is narrowband, why increasing
- the 14N RF field produces a Bloch-Siegert shift on the observed 1H
- resonance, and why B1 inhomogeneity turns that shift into broadening.
- The "offset-tolerant low-power" trace represents the design goal of
- Bloch-Siegert-aware robust optimal control or low-power amplitude-
- modulated decoupling. Replace the effective coefficients below with a
- more detailed Hamiltonian model for quantitative work.
- Magnetic field corresponding to 800 MHz 1H
- Three inequivalent 14N sites around the decoupler carrier

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `gaussian_weights()`, `decoupling_profile()`, `proton_line()`, `subplot()`, `xline()`, `site_offsets_hz()`, `num2str()`, `normalise_line()`, `site_weights()`, `b1_weights()`, `b1_scales()`, `lorentzian()`, `linewidth()`, `freq_axis_hz()`, `idx()`.
