# examples/imaging/gradient_image_1d.m

- Signature: `gradient_image_1d()`

## Purpose

Simulate a 1D imaging experiment with a hard pulse in the presence of diffusion and flow. Calculation time: seconds. Ahmed Allami and Ilya Kuprov.

## Physical / mathematical content

- The spin system contains one `1H` spin at 5.9 T with a scalar chemical shift of 1.0. It uses `t1_t2` relaxation with diagonal terms retained, zero equilibrium, and R1 and R2 rates of 0.5 and 2.0.
- The 0.30-long spatial domain has 100 points and a third-order periodic derivative. Uniform initial `Lz` and detection `L+` phantoms are specified, along with a zero relaxation phantom. Flow is `1e-2` at every spatial point, and diffusion is `5e-6`.

## Numerical / algorithmic content

- The simulation uses the `sphten-liouv` formalism without basis approximation and runs the `basic_1d_hard` sequence through `imaging` with a readout gradient amplitude of `30e-3`.
- Acquisition settings specify a sweep of 500000, 128 points, zero filling to 128 points, a zero offset, a `kHz` axis, and axis inversion. The FID receives `sqsin` apodisation; the plotted image is the real part of its shifted Fourier transform.

## Implementation structure

- Create the spin system and basis, define sequence and spatial parameters, run `imaging(spin_system,@basic_1d_hard,parameters)`, apodise the FID, Fourier-transform it, and display the result with `plot_1d`.
