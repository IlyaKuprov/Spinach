# examples/imaging/spiral_2d.m

- Signature: `spiral_2d()`
- Source: [examples/imaging/spiral_2d.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/spiral_2d.m)

## Experiment

This driver sets up two-dimensional spiral k-space imaging. Its comment estimates calculation time in minutes; that is not a measured run time here.

## Model and acquisition parameters

The spin system is one 1H with magnetic-induction parameter 5.9 and zero chemical shift; the source does not label units for these field/shift values. It uses diagonal T1/T2 relaxation, zero equilibrium, and R1=30 and R2=70; the source supplies no units for these rates. It disables the algorithmic option named 'pt' and uses the sphten-liouv formalism without a basis approximation.

The field of view is entered as dims=[0.30 0.25] and the grid as npts=[108 90]; the source does not annotate length units for dims. The acquisition settings are grad_amp=5e-2, spiral_dur=5e-3, spiral_npts=3000, spiral_frq=5e4, and t_echo=0.050. The driver does not explicitly label units for these parameter values.

R1 and R2 operators are obtained with rlx_t1_t2(spin_system), and the corresponding phantom maps R1Ph and R2Ph are loaded from [letter_a.mat](https://github.com/IlyaKuprov/Spinach/blob/main/etc/phantoms/letter_a.mat). Initial-state and receiver-state operators are Lz and L+, with uniform phantom entries. Both flow components are zero and diffusion is zero.

## Computation and observables

The driver calls imaging(spin_system,@spiral,parameters). It then derives readout and phase-encode gradient timing values from spiral_dur and spiral_npts, assigns the readout and phase-encode gradient amplitudes from grad_amp, and displays three panels: the recorded image, the R1 phantom, and the R2 phantom. The source contains no numerical image metrics or saved output values; no simulation was run for this page.

## Scope

This is a simulated two-dimensional acquisition using supplied relaxation maps and zero flow/diffusion settings. The source does not report a comparison against a measured image or establish reconstruction accuracy.

## Attribution

Ahmed Allami; Ilya Kuprov.
