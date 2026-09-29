# experiments/nmr_liquids/hmbc.m

**Canonical MATLAB source:** [experiments/nmr_liquids/hmbc.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/hmbc.m) · [Spin Dynamics Wiki: hmbc.m](https://spindynamics.org/wiki/index.php?title=hmbc.m)

## Purpose

Magnitude-mode HMBC acquisition. The implementation combines the supplied Hamiltonian, relaxation, and kinetics matrices into its Liouvillian. For natural-abundance HMBC, use Spinach isotope dilution via [dilute.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/dilute.m) to represent the isotopomer mixture, as recommended in the source header.

## Input contract

Call `fid=hmbc(spin_system,parameters,H,R,K)`. `parameters.sweep=[F1 F2]` contains the sweep widths in Hz; `parameters.npoints=[F1 F2]` contains the time-domain point counts. `parameters.spins={F1,F2}` gives nuclei in dimension order (source example: `{'15N','1H'}`). `parameters.J` is the primary scalar coupling in Hz. `parameters.delta_b` is the user-supplied delay described as `delta_2` from the cited paper; the source notes the authors' recommendation of 60e-3 s, not a built-in default. The source requires `sphten-liouv` and numeric, same-sized matrix inputs `H`, `R`, and `K`.

## Sequence and acquisition

The initial state and detection operator are `Lz` and `L+` on F2. After a +90° F2 x-pulse, the sequence evolves for `delta_a=abs(1/(2*J))`, applies a +90° F1 x-pulse, then evolves for `parameters.delta_b`. The next F1 pulse is the difference of +90° and −90° branches, followed by explicit coherence selection of +1 on F1. During sampled F1 evolution, the sequence splits the evolution around a π F2 x-pulse; it then applies a +90° F1 x-pulse, evolves for `delta_a` again, and observes F2. The `delta_a` delays are in seconds for `J` in Hz; `delta_b` is also supplied in seconds. This records the source's selection and pulse operations without assigning a longer-range transfer mechanism beyond the code.

## Output

`fid` is a two-dimensional magnitude-mode FID with `npoints(1)` F1 evolution samples and `npoints(2)` F2 detection samples.

## References

- [HMBC sequence reference](https://doi.org/10.1021/ja00268a061)
- [HMBC sequence reference](https://doi.org/10.1016/0022-2364(88)90172-2)
