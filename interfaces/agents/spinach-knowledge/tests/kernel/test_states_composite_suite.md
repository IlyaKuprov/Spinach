# tests/kernel/test_states_composite_suite.m

- Signature: `result=test_states_composite_suite()`

## Purpose

Tests composite state generators in `kernel/states`. The state helper functions must produce the expected density operators and partner-state descriptors.

## Physical / mathematical content

- Constructs two-spin singlet and triplet projectors from wavefunctions in the Zeeman product basis: `alpha`, `beta`, `sing=(ab-ba)/sqrt(2)`, and `trip_zero=(ab+ba)/sqrt(2)`.
- Checks that `four_spin_states(...,'S(x)S')` gives the tensor product of singlet projectors on spins 1–2 and 3–4.

## Numerical / algorithmic content

- Uses one-, two-, four-, and three-proton test spin systems in the `zeeman-hilb` formalism; the one-proton unit-state test also uses `zeeman-liouv`.
- Compares generated states with explicit reference matrices or vectors using `test_close`. Checks partner-state descriptor order with `test_true`.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

1. Check that the one-proton Hilbert-space `unit_state` is the `2`-by-`2` identity and that its Zeeman-Liouville counterpart is the normalised vectorised identity.
2. Check the two-proton singlet, triplet-up, triplet-zero, and triplet-down projectors against projectors built from product-basis wavefunctions.
3. Check the four-proton `S(x)S` state against `kron(singlet_pair,singlet_pair)`.
4. With spin 2 fixed at `L+`, enumerate spins 1 and 3 over `E` and `Lz`. Check the descriptor order and compare each partner-state element with a direct `state()` call.