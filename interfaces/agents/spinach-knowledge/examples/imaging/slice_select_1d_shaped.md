# examples/imaging/slice_select_1d_shaped.m

- Signature: `slice_select_1d_shaped()`

## Purpose

Shows 1D slice selection with a Gaussian-shaped RF pulse while diffusion and flow are present. The source estimates seconds of runtime.

## Model and sequence

The model is one `1H` at 5.9 T with zero chemical shift, diagonal T1/T2 relaxation, zero equilibrium, and rates `r1 = 30`, `r2 = 70`. It uses the `sphten-liouv` formalism without a basis approximation. The 0.30 m sample has 500 points; slice-selection and readout gradient amplitudes are each `30 mT/m`.

The 50-step RF pulse has total duration `0.5e-4 s`, frequency `+100 kHz`, phase `pi/2`, and a Gaussian amplitude profile scaled by `2*pi*20000`. The source sets the flow field to `1e-2` and diffusion to `5e-6`.

## Computation and output

The pulse sequence is evaluated with `imaging(spin_system,@slice_select_1d,parameters)`. The signal receives square-sine apodisation and a real, shifted Fourier transform before plotting as a 1D profile.

## Attribution

Ahmed Allami; Ilya Kuprov.
