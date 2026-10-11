# experiments/nmr_solids/wise.m

MATLAB source: [experiments/nmr_solids/wise.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_solids/wise.m)

WISE (WIdeline SEparation) is described in the source as a powder-MAS heteronuclear correlation experiment. For the common `1H`-`13C` case, it links proton line shapes in one dimension with carbon chemical shifts in the other. Reference: [10.1021/ma00038a037](https://doi.org/10.1021/ma00038a037).

## Inputs and source-defined sequence

`fid=wise(spin_system,parameters,H,R,K)` expects the working spins in order: high-gamma channel first and low-gamma channel second (the source example is `{'1H','13C'}`). It receives `H`, `R`, and `K`, combines them as `L=H+1i*R+1i*K`, and uses `spc_dim` from the context to extend the spin-control operators. The source's preparation step analytically decouples the second spin from `rho0` before building the pulse sequence.

Required fields are `spins`, `hi_pwr`, `cp_pwr`, `cp_dur`, `rho0`, `coil`, `sweep`, `npoints`, and `spc_dim`. `hi_pwr` is the high-gamma-channel RF amplitude in Hz. `cp_pwr` supplies two positive amplitudes in Hz, one per channel; `cp_dur` is the contact time in seconds. `sweep` and `npoints` each have two entries for F1 and F2; sweep widths are positive and in Hz, with dwell times `1./sweep`.

The source applies high-power 90-degree pulses along X and Y on the first channel, using `1/(4*hi_pwr)` seconds. It evolves both quadratures for `npoints(1)-1` F1 steps. The contact generator adds `-2*pi*cp_pwr(1)*Hy` on channel 1 and `+2*pi*cp_pwr(2)*Cx` on channel 2; contact evolution is for `cp_dur`. It then analytically decouples the first channel during acquisition and detects the F2 evolution on `coil` for `npoints(2)-1` steps.

## Output and limits

The return value is a structure with `fid.cos` and `fid.sin`, the cosine and sine States-quadrature components. Each contains F2 observable samples for the F1 trajectory states: `npoints(2)` by `npoints(1)`, including the initial sample on each dimension. The source does not construct a MAS rotor or powder-orientation sweep; those dynamics are represented by the context-supplied matrices and spatial dimension. The implementation describes the transfer/control sequence, not a processed spectrum.

https://spindynamics.org/wiki/index.php?title=wise.m
