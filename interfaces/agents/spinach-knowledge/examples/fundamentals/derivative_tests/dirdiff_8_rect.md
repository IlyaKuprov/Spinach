# examples/fundamentals/derivative_tests/dirdiff_8_rect.m

- MATLAB implementation: [examples/fundamentals/derivative_tests/dirdiff_8_rect.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/derivative_tests/dirdiff_8_rect.m)

- Signature: `dirdiff_8_rect()`

## Question and callable context

This no-argument example checks whether selected columns of the analytical phase Hessian returned by `grape_phase` agree with centred finite differences of its returned gradient for GRAPE with the rectangle integrator. Run `dirdiff_8_rect()` in a Spinach MATLAB environment with the mapped helper `dirdiff_test_system.m` available. The source contains error paths for failed comparisons; this page does not claim that the script was executed or that any comparison passed.

## Model and control setup

The script repeats the check for `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb`, obtaining the drift Hamiltonian `H` and the 13C states and operators from `dirdiff_test_system`. It assigns `H` as the drift, maps all four controls to channel 1, and uses control operators `Lx`, `Ly`, `0.4*Lx+0.2*Ly`, and `0.7*Ly-0.1*Lx`. The initial states are `Sx`, `Sy`, and `Sz`; their targets are `-Sz`, `Sy`, and `Sx`.

The configured power levels are `2*pi*linspace(50e3,70e3,10)`. The source sets `method='newton'`, `max_iter=1000`, no plotting options, and `integrator='rectangle'`; the demonstrated calls request values and derivatives directly from `grape_phase`. There are five intervals with `pulse_dt=12.8e-6*ones(1,5)` and amplitude rows `ones(1,5)` and `0.8+0.1*(1:5)`.

## Numerical comparison and scope

For each formalism the source draws one unseeded `2`-by-`5` phase array as `randn(2,5)/3` and uses `h=1e-5`. It obtains the analytical Hessian, then perturbs linear waveform entries 1, `end`, and 5 in turn. For each entry it estimates the corresponding Hessian column from the centred gradient difference `hess_num=(grad_forw-grad_back)/(2*h)`. The source accepts a column only if `norm(hess_anl(:,i)-hess_num,1)<1e-5*norm(hess_num,1)`; otherwise it raises an error identifying the formalism and column. Its middle-column comment refers to the specific linear index 5.

This samples three of the ten waveform entries for one random starting array per formalism; it is not a full-Hessian comparison or a sweep over random seeds and model settings. No observed pass/fail output is recorded here. Source: [examples/fundamentals/derivative_tests/dirdiff_8_rect.m](../../../../../../examples/fundamentals/derivative_tests/dirdiff_8_rect.m).