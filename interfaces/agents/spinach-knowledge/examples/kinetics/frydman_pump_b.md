# examples/kinetics/frydman_pump_b.m

- Signature: `frydman_pump_b()`

## Purpose

This is the source's Figure 9 water-exchange spin-lock pump example, attributed to Lucio Frydman. The source cites [the JMR article](https://doi.org/10.1016/j.jmr.2021.107083) and estimates a calculation time of seconds. It simulates a sequence and plots trajectories; it does not itself report a fitted exchange rate or a numerical endpoint.

## Spin and exchange model

The model has four peptide spins in the order H, N, C-alpha, and carbonyl C, plus 100 water protons. It sets `sys.magnet=11.7` (the source gives no unit beside this field parameter), all water-proton shifts to zero, and the listed H-N and N-carbonyl-C scalar couplings to -45 and 8. The source does not annotate units for these coupling entries. It assigns both T1 and T2 values in seconds: H 0.2722, N 0.8, C-alpha 2, carbonyl C 2, and water 0.2994. Relaxation is `t1_t2`, diagonal terms are retained, equilibrium is `IME`, and temperature is 298.

The intermolecular flux-rate matrix couples the amide N (spin 2) to the first water proton (spin 5) with a value of 1000 in each direction. The water-water block is filled with 1e4. These inputs are named exchange rates but have no units annotated in the source; the second assignment is a full block, not a reported measurement of a single water-water rate.

## Sequence and propagation

The function uses the `sphten-liouv` formalism with `IK-1`, scalar-coupling connectivity, `inter_level=4`, and `prox_level=1`. It starts from isotropic thermal equilibrium and executes ten loops. Each loop applies the source-described flip-lock-backflip-lock-backflip pump: the sequence performs forward 90-degree flips, propagates through two cross-polarisation periods, applies the listed proton, carbon, and nitrogen back-flips, then destroys residual coherences with `homospoil`. The two effective DIPSI Hamiltonians are built with different Zeeman and coupling settings for the two periods.

The raw duration and sampling inputs are `cp_dur=[11e-3 53e-3]`, `cp_npt=[11 53]`, and `nloops=10`. The output trajectory samples both periods in each loop. The plot labels its horizontal axis as time in milliseconds.

## Observable and scope

The trajectory is projected onto H, N, and carbonyl-C longitudinal and transverse operators (`Lz` and `Lx`). Three stacked plots show the paired Z/X expectation-value traces against time.

This is a forward sequence simulation with source-specified relaxation and exchange inputs. The page does not infer a fitted rate or quote a trajectory value not printed by the source. See the [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/frydman_pump_b.m) for the complete setup and pulse implementation.
