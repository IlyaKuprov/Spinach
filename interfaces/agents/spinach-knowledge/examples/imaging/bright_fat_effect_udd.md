# examples/imaging/bright_fat_effect_udd.m

- MATLAB implementation: [examples/imaging/bright_fat_effect_udd.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/bright_fat_effect_udd.m)

## Purpose

This example models the bright-fat effect under a UDD echo train. Its source comment describes greater magnetisation losses in MRI experiments on J-coupled systems as coherences are lost in the depths of Hilbert space. The source estimates minutes of simulation time and notes it is faster with a Tesla V100 GPU.

## Spin system and basis

The model is the six-spin, two-molecule setup also used in the CPMG example: six 1H spins at `sys.magnet=3.0`, with spins 1-3 assigned to molecule A and 4-6 to molecule B. Each molecule has chemical-shift entries `{1.0,2.0,3.0}`; the source states no units for these entries or the magnetic-induction setting. Couplings within A are zero; within B, pairs (4,5), (4,6) and (5,6) receive scalar values 11, 17 and 23, with no unit stated. The two chemical parts have an empty reaction list (no chemical exchange) and concentrations `[1,1]`.

The basis is `sphten-liouv` with no approximation; path tracing is disabled. GPU enablement is commented out, and relaxation phantom and operator lists are empty.

## Sequence and imaging setup

The sequence uses the 1H channel, no decoupling, zero offset, 48 pulses, and `dec_time=80e-3`; the source does not give the time unit. Geometry is `dims=[0.30,0.25]` with `npts=[100,200]` and derivative setting `{'period',3}`. It loads the left and right bright-fat phantoms from `../../etc/phantoms/bright_fat_left.mat` and `../../etc/phantoms/bright_fat_right.mat`, assigns `{1-left,1-right}` as phantom initial states, uses Lz states for the two molecule spin groups, and detects the 1H Lx state with a uniform coil phantom. Flow fields are zero and `diff=0`.

The distinctive sequence call is `imaging(spin_system,@udd_dec,parameters)`. The script plots `surf(abs(mri))` with the X direction reversed and pixel-labelled axes, under the title 'Bright fat effect under UDD echo train'. It defines no numerical image result in the page.

The detection vector uses unweighted `coil_state`; initial magnetisation uses concentration-weighted `state`.
