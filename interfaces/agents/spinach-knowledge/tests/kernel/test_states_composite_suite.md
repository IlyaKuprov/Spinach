# tests/kernel/test_states_composite_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_states_composite_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_states_composite_suite.m)

## Purpose

Regression test suite for the composite state generators in `kernel/states`. The suite verifies that state helper functions produce the textbook density operators, covering:

- Unit-state normalisation in both Zeeman-Hilbert and Zeeman-Liouville formalisms.
- Two-spin singlet and triplet projectors.
- Four-spin product states.
- Partner-state enumeration.

## Behaviour

The function announces the test target with `fprintf('TESTING: Composite state generators\n')` and initialises a regression test result via `new_test_result('kernel/states_composite_suite', 'Composite state generators', 'state helper functions must produce the textbook density operators.')`.

The suite then performs the following checks:

1. **Unit state, Zeeman-Hilbert:** Builds a one-proton spin system (`sys.isotopes={'1H'}`, `bas.formalism='zeeman-hilb'`, `bas.approximation='none'`) using `test_spin_system`, then compares `unit_state(spin_system)` against `speye(2)` with tolerances `1e-15` (absolute and relative), asserting that the Hilbert-space unit state is the identity density matrix.

2. **Unit state, Zeeman-Liouville:** Rebuilds the spin system with `bas.formalism='zeeman-liouv'`, forms the vectorised normalised identity `unit=speye(2); unit=unit(:); unit=unit/norm(unit,2)`, and compares `unit_state(spin_system)` against this vector with tolerances `1e-15`, asserting that the Zeeman-Liouville unit state is the normalised vectorised identity.

3. **Singlet and triplet projectors:** Builds a two-proton spin system (`sys.isotopes={'1H','1H'}`, Zeeman-Hilbert formalism). Using the Zeeman product basis vectors `alpha=[1;0]`, `beta=[0;1]`, and their Kronecker products `aa`, `ab`, `ba`, `bb`, the suite constructs `sing=(ab-ba)/sqrt(2)` and `trip_zero=(ab+ba)/sqrt(2)` and checks:
   - `singlet(spin_system,1,2)` against `sing*sing'` with tolerances `1e-14`.
   - `[TU,T0,TD]=triplet(spin_system,1,2)` against `aa*aa'`, `trip_zero*trip_zero'`, and `bb*bb'` respectively, each with tolerances `1e-14`.

4. **Four-spin product states:** Builds a four-proton spin system (`sys.isotopes={'1H','1H','1H','1H'}`, Zeeman-Hilbert formalism) and checks `four_spin_states(spin_system,[1 2 3 4],'S(x)S')` against `kron(singlet_pair,singlet_pair)` where `singlet_pair=sing*sing'`, with tolerances `1e-13`, asserting that `S(x)S` is the tensor product of singlets on spins 1-2 and 3-4.

5. **Partner-state enumeration:** Builds a three-proton spin system (`sys.isotopes={'1H','1H','1H'}`, Zeeman-Hilbert formalism). It calls `[A,descr]=partner_state(spin_system,{{'L+',2}},{{{'E','Lz'},[1 3]}})`, enumerating two partner spins (1 and 3) that may each be `E` or `Lz` while spin 2 is `L+`. The expected descriptor order is:

   ```matlab
   {{'E','L+','E'};
    {'Lz','L+','E'};
    {'E','L+','Lz'};
    {'Lz','L+','Lz'}}
   ```

   The suite checks via `test_true` that `descr` matches this expected order, asserting that `partner_state()` enumerates the Cartesian product of allowed partner states. Each enumerated state `A{n}` is then compared against a direct `state(spin_system,expected_descr{n},{1,2,3})` call with tolerances `1e-14`, asserting that each partner-state descriptor maps to the corresponding direct product state.

All comparisons use `test_close` (or `test_true` for the descriptor order check) to accumulate results and explanatory messages into the returned regression test result.

## Inputs and outputs

```matlab
result = test_states_composite_suite()
```

**Outputs:**

- `result` — regression test result with explanatory messages.

**Inputs:**

- None.

## References

- [Spinach GitHub repository — test_states_composite_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_states_composite_suite.m)
