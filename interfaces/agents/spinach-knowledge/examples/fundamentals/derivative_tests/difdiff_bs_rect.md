# examples/fundamentals/derivative_tests/difdiff_bs_rect.m

- MATLAB implementation: [examples/fundamentals/derivative_tests/difdiff_bs_rect.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/derivative_tests/difdiff_bs_rect.m)

- Signature: `difdiff_bs_rect()`

## Purpose

Checks selected Cartesian GRAPE fidelity-gradient entries against centred finite differences when the rectangular integrator and Bloch–Siegert corrections are enabled. It checks the first waveform entry, last entry, and one midpoint entry; it is not an element-by-element gradient comparison.

## Spin system and control setup

The four-spin system sets `sys.magnet=10.2`, isotopes `{'1H','1H','13C','13C'}`, scalar Zeeman values `{1.5,2.0,30.0,40.0}`, and scalar couplings 1–2: `7.0`, 1–3: `150`, 2–4: `150`, and 3–4: `50`. It uses the `sphten-liouv` basis with `approximation='none'`. The normalised initial state is the singlet on spins 1 and 2; the normalised target is `Lz` on spin 4.

The controls are `Lx` and `Ly` for each of `1H` and `13C`, with channels `[1,1,2,2]`. The drift is the NMR Hamiltonian. The control configuration sets offsets `{1050,5285}` with corresponding `Lz` offset operators, `pwr_levels=2*pi*500`, `integrator='rectangle'`, 100 pulse intervals of `1.5e-4`, and `max_iter=1000`. Bloch–Siegert correction is enabled for isotopes `{'1H','13C'}`. The source does not state units for these numeric settings.

## Gradient comparison

After `optimcon`, the script draws one random `4×100` guess waveform as `randn(...)/10` and uses finite-difference increment `h=1e-5`. It obtains the analytical gradient from `grape_xy`, then perturbs one waveform entry at a time by `±h`; the numerical derivative is `(fid_forw(1)-fid_back(1))/(2*h)`. The sampled entries are the first and last waveform entries and the entry at row `ceil(4/2)=2`, column `ceil(100/2)=50`.

At each entry the relative error is `abs(grad_anl-grad_num)/max(abs(grad_num),eps('double'))`. An entry passes only when that value is below `5e-6`. Individual failures include the analytical gradient, finite-difference gradient, and relative error; after all three checks, any failure causes the script to raise an error. Although `max_iter` is configured, this script's shown comparison calls `grape_xy` and does not run an optimisation loop.

## Scope

The waveform is randomised without a source-level seed, and only three of its 400 entries are checked in a run. These checks do not establish the accuracy of every gradient component or report an optimiser outcome. No units are added to the source's magnetic-field, offset, power, or pulse-interval values.
