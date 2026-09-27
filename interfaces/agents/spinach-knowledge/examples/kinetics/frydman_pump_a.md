# examples/kinetics/frydman_pump_a.m

- Signature: `frydman_pump_a()`

## Purpose

Lucio Frydman's water-exchange spin-lock pump, Figure 2 from https://doi.org/10.1016/j.jmr.2021.107083. The script calculates the H/N trajectory for a single flip-and-lock cycle; the source estimates seconds.

## Physical / mathematical content

The model contains a four-spin peptide system (H, N, Cα, C′) and 20 water protons. An intermolecular flux matrix connects the amide proton to the first water proton at 10 Hz and couples the water pool at 10⁴ Hz. The peptide H–N scalar coupling is −45 Hz. Relaxation uses the source's diagonal T1/T2 rates and isotropic equilibrium at 298 K.

## Numerical / algorithmic content

The callback removes the Zeeman terms and treats the H–N, N–C, and H–C couplings as strong for the effective spin-lock Hamiltonian. After a nitrogen crusher, 90° H and N flips prepare the state; the 100 ms CP period is propagated with the Krylov trajectory method using 100 points.

## Implementation structure

- Sets B₀ = 11.7 T, 20 water protons, and the peptide isotope/coupling model.
- Uses `t1_t2` relaxation with diagonal retention and `IK-1` basis settings.
- Sets NH–water and water-pool exchange rates to 10 and 10⁴ Hz, respectively.
- Calls `liquid` with `frydman_pump`; plots H and N Z/X expectation values versus time.
