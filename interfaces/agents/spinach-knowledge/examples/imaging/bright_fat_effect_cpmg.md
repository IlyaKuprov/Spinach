# examples/imaging/bright_fat_effect_cpmg.m

- MATLAB implementation: [examples/imaging/bright_fat_effect_cpmg.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/bright_fat_effect_cpmg.m)

## Purpose

This example models the bright-fat effect under a CPMG echo train. Its source comment describes greater magnetisation losses in MRI experiments on J-coupled systems as coherences are lost in the depths of Hilbert space. The source estimates minutes of simulation time and notes it is faster with a Tesla V100 GPU.

## Spin system and basis

The system has six 1H spins at `sys.magnet=3.0`. Spins 1-3 are molecule A and spins 4-6 are molecule B. Each molecule has chemical-shift entries `{1.0,2.0,3.0}`; no units are stated for these entries or for the magnetic-induction setting. Scalar couplings within A are set to zero; within B the source assigns 11, 17 and 23 to pairs (4,5), (4,6) and (5,6), without stating units. The two chemical parts have an empty reaction list (no chemical exchange) and concentrations `[1,1]`.

The basis uses `sphten-liouv` with no approximation, and path tracing is disabled. GPU enablement is present only as a commented-out line. Relaxation phantom and operator lists are empty.

## Sequence and imaging setup

The sequence uses `parameters.spins={'1H'}`, no decoupling, zero offset, 48 pulses and `dec_time=80e-3`. The source does not specify a unit for `dec_time`. Sample geometry is `dims=[0.30,0.25]` on `npts=[100,200]`, with derivative setting `{'period',3}`. It loads `left` and `right` from `../../etc/phantoms/bright_fat_left.mat` and `../../etc/phantoms/bright_fat_right.mat`, uses phantom initial states `{1-left,1-right}` and Lz states for the two three-spin groups, and detects the 1H Lx state with a uniform coil phantom. Flow fields `u` and `v` are zero and `diff=0`.

The simulation call is `imaging(spin_system,@cpmg_dec,parameters)`. The source plots `surf(abs(mri))`, reverses the X direction, and labels the axes in pixels. No numerical image result is stated in the source page.

The detection vector uses unweighted `coil_state`; initial magnetisation uses concentration-weighted `state`.
