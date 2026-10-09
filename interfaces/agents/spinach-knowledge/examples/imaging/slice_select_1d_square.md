# examples/imaging/slice_select_1d_square.m

- Signature: `slice_select_1d_square()`
- Source: [examples/imaging/slice_select_1d_square.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/slice_select_1d_square.m)

## Experiment

This driver runs one-dimensional slice selection with a rectangular RF pulse while the sample has both flow and diffusion. Its source comment estimates a runtime of seconds; that is a source estimate, not a measured timing here.

## Model and encoded parameters

The spin system contains one 1H with magnetic-induction parameter 5.9 and zero chemical shift; the source does not label a unit for the field value. It uses the sphten-liouv formalism with no basis approximation. The T1/T2 relaxation model keeps diagonal terms, sets the equilibrium state to zero, and assigns R1=30 and R2=70; the source gives no units for these rate values.

The sample geometry is entered as dims=0.30 with 100 points; the source does not annotate a unit for dims. The spatial derivative option is {'period',3}. Readout settings are sweep=500000, 128 points, axis_units='kHz', and invert_axis=1. Slice-selection and readout gradient amplitudes are each 30e-3, with no units stated in the driver. Diffusion is 5e-6 and the one-dimensional flow field is 1e-2 at each sample point; their units are likewise not stated.

The rectangular RF pulse is represented by 50 segments over a total pulse_time of 0.5e-4. The source sets the pulse frequency to +100e3, pulse power to `2*pi*3000`, and RF phase to pi/2, without annotating units for frequency or power. It supplies segment-wise frequency, amplitude, and duration lists, and sets max_rank=3.

The initial-state and coil phantoms are uniform, with Lz preparation and L+ detection operators. The relaxation phantom is zero and uses relaxation(spin_system) as its operator.

## Computation and observable

After creating the spin system and basis, the driver calls imaging(spin_system,@slice_select_1d,parameters). It applies square-sine apodisation, takes a shifted Fourier transform, keeps its real part, and plots the resulting one-dimensional profile. The source specifies the calculation but does not give numerical profile values or a saved result; no simulation was run for this page.

## Scope

This example demonstrates a one-dimensional pulse-and-gradient acquisition with flow and diffusion enabled. The source does not report an experimental comparison or quantify slice-profile performance.

## Attribution

Ahmed Allami; Ilya Kuprov.
