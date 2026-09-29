# kernel/utilities/merge_inp.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/merge_inp.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/merge_inp.m)

## Purpose

Merges multiple `sys` and `inter` structures into one. Useful for setting up chemical kinetics simulations where the molecules come from different DFT calculations.

## Behaviour

- Syntax: `[sys,inter]=merge_inp(sys_parts,inter_parts)`.
- Consistency is enforced first: both inputs must be row cell arrays of structures, `sys_parts` and `inter_parts` must have the same number of elements, at least one subsystem must be supplied, and every `sys` structure must contain an `isotopes` subfield.
- Spin counts per subsystem are obtained from `numel(x.isotopes)` on each `sys` entry; these sizes drive the spin-index offsets used when merging index-based fields.
- Extensive fields are concatenated with spin and chemical subsystem indices offset as appropriate; non-extensive fields (magnet, temperature, relaxation settings, etc.) must be identical in all subsystems that supply them, and any difference is treated as an error.
- `sys` fields merged as common values that must be equal: `magnet`, `output`, `scratch`, `enable`, `disable`, `tols`, `parallel`, `parprops`.
- `sys` fields merged as row cell arrays: `isotopes`, `labels`.
- `inter` fields merged as common values that must be equal: `relaxation`, `rlx_keep`, `rlx_dfs`, `equilibrium`, `temperature`, `damp_rate`, `srfk_tau_c`, `nz_shift`, `nz_onshell`, `weiz_r1e`, `weiz_r2e`, `weiz_r1n`, `weiz_r2n`, `nott_r1e`, `nott_r2e`, `nott_r1n`, `nott_r2n`, `pbc`.
- `tau_c` and `order_matrix` are merged per chemical species when a `chem` group is present and every subsystem's `chem` contains a `parts` subfield; otherwise they are merged as common values that must be equal.
- `coordinates` is merged as a column cell array (vertical concatenation).
- Per-spin relaxation rate fields `r1_rates`, `r2_rates`, `lind_r1_rates`, `lind_r2_rates` are vertically concatenated as column vectors.
- Square arrays over spin pairs (`srfk_mdepth`, `weiz_r1d`, `weiz_r2d`) are merged block-diagonally for numeric arrays; cell-array forms are merged into a larger cell array with each subsystem's block placed on the diagonal, and mixed cell/non-cell types across subsystems are an error.
- Spin index vectors (`srsk_sources`) are reshaped to row vectors and concatenated with per-subsystem offsets (`cumsum` of preceding part sizes). Cell arrays of spin index vectors (`ignore`) are merged the same way, with each vector's indices shifted by the subsystem offset.
- Nested groups `zeeman`, `coupling`, `giant`, `suscept` and `chem` must each be present in all subsystems or in none; each is merged field by field and then removed from the subsystem copies:
  - `zeeman`: `matrix`, `scalar`, `eigs`, `euler` merged as row cell arrays.
  - `coupling`: `matrix`, `scalar`, `eigs`, `euler` merged as square (block-diagonal) arrays.
  - `giant`: `coeff`, `euler` merged as row cell arrays.
  - `suscept`: `chi`, `xyz` merged as row cell arrays.
  - `chem`: `parts` merged as cell arrays of spin index vectors with offsets; `rates` and `flux_rate` merged as square arrays; `concs` merged as a row cell array; `flux_type`, `rp_theory`, `rp_rates` merged as common values that must be equal; `rp_electrons` merged as spin index vectors with offsets.
- Any field remaining in `sys_parts` or `inter_parts` after all merges triggers an error (`unhandled subfield in sys.` / `unhandled subfield in inter.`); similarly, unhandled subfields inside a nested group trigger an error naming the group.
- Coordinates and susceptibility centres from all subsystems are assumed to refer to one common frame of reference; spin index lists are returned as row vectors.

## Inputs and outputs

**Inputs**

- `sys_parts` — a row cell array of `sys` structures to be merged.
- `inter_parts` — a row cell array of `inter` structures to be merged.

**Outputs**

- `sys` — resulting `sys` structure.
- `inter` — resulting `inter` structure.

## References

- Spinach Wiki: [merge_inp.m](https://spindynamics.org/wiki/index.php?title=merge_inp.m)
- Source code contact: ilya.kuprov@weizmann.ac.il, a.acharya@soton.ac.uk
