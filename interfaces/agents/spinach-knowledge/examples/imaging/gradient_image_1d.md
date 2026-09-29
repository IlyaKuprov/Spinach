# examples/imaging/gradient_image_1d.m

A one-dimensional hard-pulse imaging example with diffusion and flow; its header gives a seconds-scale calculation time and credits Ahmed Allami and Ilya Kuprov. [Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/gradient_image_1d.m)

## Model and sequence

The model is one `1H` spin with `sys.magnet=5.9` and scalar chemical shift `1.0`; the example does not annotate units for these values. Relaxation is `t1_t2` with diagonal terms retained, zero equilibrium, `R1=0.5` and `R2=2.0`. It uses `sphten-liouv` without basis approximation. The sequence call is `imaging(spin_system,@basic_1d_hard,parameters)`, with no decoupled spins, offset `0`, sweep value `500000`, 128 acquired points and 128-point zero-fill; the displayed axis is configured in kHz and inverted. A readout-gradient amplitude of `30e-3` is supplied without a unit comment.

The sample uses `dims=0.30`, 100 points and `{'period',3}` differentiation. The relaxation phantom is all zeros, with the relaxation operator from `relaxation(spin_system)`. Initial and coil spatial phantoms are uniform ones; their states are `Lz` and `L+`, respectively. Flow is `u=1e-2` at every point and diffusion is `5e-6`; the example does not annotate units for these values or the sample length.

## Output and interpretation

The imaging result is square-sine apodised, then transformed with a centred 1D FFT; the real spectrum is passed to `plot_1d`. The source gives the display axis unit (`kHz`), but does not label the gradient amplitude, flow, diffusion, or geometry units. It specifies the hard-pulse sequence helper but no independent RF amplitude or duration.
