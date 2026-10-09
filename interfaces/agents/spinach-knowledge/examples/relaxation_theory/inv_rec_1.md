# examples/relaxation_theory/inv_rec_1.m

- MATLAB implementation: [examples/relaxation_theory/inv_rec_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/inv_rec_1.m)

Source: [examples/relaxation_theory/inv_rec_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/inv_rec_1.m)

## Purpose

Simulates one proton's longitudinal signal after an inversion pulse, then plots that signal during a one-second recovery. The source describes it as a simple inversion-recovery example and estimates seconds of calculation time.

## Spin and relaxation model

The model has one 1H spin at 14.1 T and sets inter.zeeman.scalar={1.5}; the source gives no unit for this scalar value. It uses the t1_t2 relaxation model with r1_rates={5.0} and r2_rates={5.0}, dibari equilibrium, secular retention (rlx_keep='secular'), and temperature value 298 (no temperature unit is stated in the source). Both configured relaxation rates are source inputs, not fitted or reported experimental measurements. The basis is the complete sphten-liouv basis.

## Pulse, detection, and plotted observable

The initial density operator is equilibrium(spin_system). The script defines an Lz detection state and an Lx pulse operator, applies a pi pulse, and evolves under L+1i*R with a 1 ms step for 1000 steps. The resulting real observable is plotted against linspace(0,1,1001) seconds with the axis labelled as the S_Z expectation value. Thus the figure is a simulated recovery trace from the configured model; no experimental curve or comparison is supplied by this source.
