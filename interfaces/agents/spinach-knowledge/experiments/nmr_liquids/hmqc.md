# experiments/nmr_liquids/hmqc.m

**Canonical MATLAB source:** [experiments/nmr_liquids/hmqc.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/hmqc.m) · [Spin Dynamics Wiki: hmqc.m](https://spindynamics.org/wiki/index.php?title=hmqc.m)

## Purpose

Magnitude-mode HMQC acquisition with indirect F1 evolution and F2 detection. The function composes the supplied Hamiltonian, relaxation, and kinetics matrices into a Liouvillian.

## Input contract

Call `fid=hmqc(spin_system,parameters,H,R,K)`. `parameters.sweep=[F1 F2]` gives sweep widths in Hz, and `parameters.npoints=[F1 F2]` gives point counts. `parameters.spins={F1,F2}` identifies nuclei in dimension order (source example: `{'15N','1H'}`). `parameters.decouple_f1` lists nuclei receiving midpoint 180° refocusing pulses during F1 (example: `{'1H'}`); `parameters.decouple_f2` lists nuclei decoupled in F2 (example: `{'15N','13C'}`). `parameters.J` is the primary scalar coupling in Hz. The source requires `sphten-liouv` and numeric, same-sized matrix inputs `H`, `R`, and `K`.

## Sequence and acquisition

The initial state is `Lz` on F2 and detection uses `L+` on F2. A +90° F2 x-pulse is followed by `delta=abs(1/(2*J))`, then a +90° F1 x-pulse and explicit selection of +1 F1 coherence. The F1 evolution is sampled in two halves, with a π pulse on each nucleus listed in `parameters.decouple_f1` between them. A further +90° F1 x-pulse precedes a second `delta` evolution; nuclei in `parameters.decouple_f2` are decoupled before F2 observation. For `J` in Hz, each `delta` is in seconds. No additional transfer-pathway interpretation is needed beyond these explicit source operations.

## Output

`fid` is a two-dimensional magnitude-mode FID with `npoints(1)` F1 evolution samples and `npoints(2)` F2 detection samples. The source recommends isotope dilution for natural-abundance experiments; see `dilute.m`.

## Reference

- [HMQC sequence reference](https://doi.org/10.1016/0022-2364(83)90241-X)
