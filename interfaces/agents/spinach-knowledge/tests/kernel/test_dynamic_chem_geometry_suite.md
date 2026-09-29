# tests/kernel/test_dynamic_chem_geometry_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_chem_geometry_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_chem_geometry_suite.m)

## Purpose

Regression test for the deterministic chemistry and geometry utility helpers in the Spinach kernel. It verifies that lattice construction, geometric measurements, coupling extraction, chemical shifts, nearest-neighbour lookup, and tensor helpers preserve exact coordinate, tensor, and metadata formulae.

## Behaviour

The test announces its target with `fprintf('TESTING: Chemistry and geometry utilities\n')` and initialises a result object via `new_test_result` for `kernel/dynamic_chem_geometry_suite`, describing the target as small chemistry and geometry helpers that must preserve exact coordinate, tensor, and metadata formulae.

A local three-spin fixture (`local_geometry_spin_system`) is built with `spin_system.sys.output='hush'` and an empty `spin_system.sys.disable` cell. It defines `nspins=3`, isotopes `{'1H','13C','15N'}`, labels `{'proton','carbon','nitrogen'}`, and `mults=[2 2 2]`. Coordinates are `{[0 0 0],[2 0 0],[0.5 0 0]}` and chemical parts are `{[1 2],3}`. The coupling matrix is a 3-by-3 cell array with `{1,2}=diag([1 2 3])` and `{2,1}=diag([4 5 6])`. Base frequencies are `2*pi*[100e6 25e6 10e6]` and offsets are `2*pi*[100 50 0]`; the Zeeman matrix cells are `basefrqs(k)*eye(3)+offsets(k)*eye(3)` for each spin. Magnetogyric ratios come from `spin('1H')`, `spin('13C')`, `spin('15N')`; `zeeman.ddscal` holds three identity tensors; `tols.hbar=1.054571628e-34` and `tols.muB=9.274009994e-24`.

The suite then performs the following checks, each recorded with `test_true` or `test_close`:

- **Cubic lattice:** calls `cubic_lattice('13C',2,2)` and asserts eight isotope entries all equal to `'13C'`; checks `inter.coordinates{2}` equals `[0 0 2]` (sub2ind ordering places the second atom one spacing along z) with tolerances `1e-15`; and checks `inter.pbc` equals `{[4 0 0],[0 4 0],[0 0 4]}` (periodic boundary vectors are spacing times the number of periods).
- **Dihedral:** `dihedral([1 0 0],[0 0 0],[0 1 0],[0 1 1])` must return `-90` degrees with tolerances `1e-12`, per the Spinach torsion-sign convention for that A-B-C-D geometry.
- **Point-cloud density:** `xyz2pd(coords,[0 1],[0 1],[0 1],2,2,2)` with `coords=[0.25 0.25 0.25;0.75 0.25 0.25;1.5 0.5 0.5]` must produce a 2-by-2-by-2 array with ones at positions `(1,1,1)` and `(2,1,1)` and zeros elsewhere; the two in-range points are counted in their x-axis bins and the out-of-range point is discarded.
- **Nearest neighbour:** `[spin_idx,dist]=nearest_spin(spin_system,1)` must return index `3` and distance `0.5` Angstrom, with tolerances `1e-15`.
- **Substance ownership:** `which_subst(spin_system,[1 2])` returns `1` (spins one and two belong to the first chemical part) and `which_subst(spin_system,3)` returns `2`.
- **Coupling extraction:** `get_coupling(spin_system,1,2)` must equal `diag([5 7 9])`, i.e. the sum of the forward and backward coupling tensor cells.
- **Chemical shifts:** `[cs_ppm,cs_hz]=chemshifts(spin_system)` must give ppm shifts `[1 2 0]` (tolerances `1e-9`/`1e-12`) and hertz shifts `[-100 -50 0]` (tolerances `1e-7`/`1e-12`); ppm shifts are one million times isotropic offsets divided by base frequencies, and hertz shifts are minus isotropic angular offsets divided by `2*pi`. `offsetof(spin_system,1)` must return `-100` with the same tolerances.
- **g-tensor:** after setting `spin_system.inter.zeeman.ddscal{1}=2*eye(3)`, `gtensorof(spin_system,1)` must equal `-ddscal{1}*gammas(1)*hbar/muB`, applying the documented `gamma*hbar/muB` scaling, with tolerances `1e-15`.
- **Isotropic shift:** `shift_iso({diag([1 2 6]),eye(3)},1,10)` must return `diag([8 9 13])` for the first tensor (replacing an isotropic part of three by ten adds seven to each diagonal component, tolerances `1e-12`) and leave the second tensor unchanged as `eye(3)` (tolerances `1e-15`).

## Inputs and outputs

Syntax:

```matlab
result = test_dynamic_chem_geometry_suite()
```

The function takes no inputs. It returns `result`, a regression test result object with explanatory messages accumulated by the `test_true` and `test_close` assertions.

## References

- Source file: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_chem_geometry_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_chem_geometry_suite.m)
