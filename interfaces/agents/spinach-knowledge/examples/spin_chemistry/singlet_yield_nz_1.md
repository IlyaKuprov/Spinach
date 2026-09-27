# examples/spin_chemistry/singlet_yield_nz_1.m

- Signature: `singlet_yield_nz_1()`

## Purpose

Compare Redfield theory with the lifetime-shifted Nakajima-Zwanzig (NZ) kernel for the magnetic-field dependence of singlet recombination yield in a triplet-born benzophenone ketyl/thiyl radical pair in a viscous ionic liquid. The source describes cage lifetimes of tens of nanoseconds and rotational correlation times of a few nanoseconds, putting `k*tau_c` near `0.1`; recombination competes with rotational decorrelation of anisotropic hyperfine couplings. The parameters represent TMPA-TFSA measurements of the Wakasa group.

## Setup

The field sweep uses 15 points from 0.5 to 50 mT, with cage drain `3e7 Hz` and correlation time `3.3 ns`. A second sweep evaluates the rising-edge yield at `3 mT` for drain rates `1e6, 3e6, 1e7, 3e7, 1e8 Hz`. Both sweeps call the cage-pair calculation with Redfield and NZ treatments.

## References

- [DOI: 10.1021/jp074331a](https://doi.org/10.1021/jp074331a)
- [DOI: 10.1039/c2cp23747d](https://doi.org/10.1039/c2cp23747d)
