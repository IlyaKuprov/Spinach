# examples/spin_chemistry/singlet_yield_nz_2.m

- Signature: `singlet_yield_nz_2()`
- Source: [examples/spin_chemistry/singlet_yield_nz_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/spin_chemistry/singlet_yield_nz_2.m)

## Purpose

Despite its function name, this example calculates the field dependence of the slowest decay rate for a micelle-confined, triplet-born benzophenone ketyl/alkyl radical pair. It compares Redfield relaxation with a lifetime-shifted Nakajima-Zwanzig (NZ) treatment. The source comments describe high-field decay as relaxation-controlled, with anisotropic hyperfine modulation draining T+/- population into the reactive S/T0 subspace.

## Spin model and relaxation

The system contains two electrons and two 1H nuclei. The scalar electron g-factors are 2.0031 and 2.0026. The electron-2/proton hyperfine tensors use principal values passed through `mt2hz([1.7 1.7 3.2])` and `mt2hz([2.1 2.1 3.3])`; their Euler angles are [0 0 0] and [pi/5 pi/3 pi/7], respectively. The Zeeman field is varied in Tesla.

Recombination is Haberkorn, for electron pair [1 2], with rates [k_rec 0]: only the singlet channel is assigned a nonzero rate, `k_rec=1e9` Hz. The model computes relaxation twice, using `redfield` and `naka-zwan`, at `tau_c=0.7e-9` s, with zero equilibrium and the lab-frame relaxation representation. For NZ, the code selects the chemical lifetime shift and sets `nz_onshell=false`. The basis is the full sphten-liouv formalism with no approximation.

## Sweep and observable

The ten field values are 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0, and 1.34 T. At each point, the script forms `L=full(H+1i*R+1i*K)`, converts its eigenvalues to decay rates with `-imag(eig(L))`, and selects the smallest rate above the 1 Hz numerical floor. It prints and plots the two slowest-mode decay-rate curves in Hz. It does not propagate an initial density operator or calculate a singlet yield; the triplet-born description is in the source comments, not an initial-state specification in the code.

The header comment says the scalar lifetime-shift parameter `k*tau_c` is near 0.35, while the explicit inputs give `1e9*0.7e-9=0.7`. The header cites [10.1246/bcsj.57.322](https://doi.org/10.1246/bcsj.57.322) and [10.1016/0009-2614(80)80325-2](https://doi.org/10.1016/0009-2614(80)80325-2) as context for representative SDS supercage parameters; these are citations present in the source comment, not independently verified experimental yields or results.
