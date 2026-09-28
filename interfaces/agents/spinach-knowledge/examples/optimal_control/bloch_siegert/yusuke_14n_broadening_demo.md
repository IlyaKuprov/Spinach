# examples/optimal_control/bloch_siegert/yusuke_14n_broadening_demo.m

- Signature: `yusuke_14n_broadening_demo()`

## Status

Historical example: the MATLAB source was removed in commit `b4f03f29`. This note describes the source at its parent revision; the example is not available at the current source path and should not be treated as a runnable current demo.

## Purpose and model

The deleted script is a qualitative reduced-model illustration, inspired by work of Nehra, Agarwal, and Nishiyama on (^{14}mathrm N) decoupling under (^{1}mathrm H) detection. It contrasts narrowband low-power CW decoupling, higher-power CW decoupling, and a hypothetical low-power offset-tolerant profile. The latter is a design target, not a pulse computed by optimal control. The source explicitly recommends replacing its effective coefficients with a more detailed Hamiltonian model for quantitative work.

At 18.8 T (approximately 800 MHz for (^{1}mathrm H)), it represents three (^{14}mathrm N) sites at −18, 0, and +18 kHz, with weights 0.30, 0.40, and 0.30, and samples a Gaussian (B_1)-scale distribution from 0.85 to 1.15. Its illustrative RF settings are 8 kHz and 20 kHz for low- and high-power CW, and 12 kHz for the offset-tolerant trace. Hand-defined response profiles determine residual broadening; a quadratic effective Bloch–Siegert shift is averaged over the (B_1) distribution and removed at its mean to show the residual line broadening.

The script plots decoupling efficiency against (^{14}mathrm N) offset and normalized predicted proton lines, and prints mean efficiencies and FWHM values. Its intended qualitative point is that a Bloch–Siegert shift becomes broadening when dispersed by (B_1) inhomogeneity. It is an analytic toy model, not a full quadrupolar/MAS calculation or a quantitative prediction. The source supplies author names but no full bibliographic citation; no additional citation details are inferred here.
