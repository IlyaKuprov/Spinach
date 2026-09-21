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
