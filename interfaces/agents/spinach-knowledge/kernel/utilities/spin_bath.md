# kernel/utilities/spin_bath.m

[MATLAB source](../../../../../kernel/utilities/spin_bath.m)

## Syntax

`[sys,inter,bath]=spin_bath(sys,inter,specification)`

## Contract

The three inputs are scalar structures. sys.magnet is authoritative in tesla; specification.magnetic is rejected. Existing isotopes, labels, coordinate rows, and native interaction-table entries are preserved when appending. Existing sites must have native coordinates and correctly sized interaction tables. An empty system needs only sys.magnet and an empty inter structure.

The explicit/lattice/random finite-bath generation formulas retain explicit coordinates and row cell vectors in Angstrom, finite repeats or bounds, occupancy, correlated isotope configurations, exclusion/selection, rigid frames, pair tensors, and candidate-site groups. Specification isotopes is a cell row; it is not silently reshaped. Random choices require an explicit local-stream seed and do not change the global random stream. No material composition, natural abundance, optical population, polarisation, or quadrupole is inferred.

tensors_Hz is a total local Zeeman tensor. Its rotated nuclear/electron values are converted to native shielding/g units; create subsequently produces the canonical rad/s tensor. Nonzero total local tensors at zero field are refused. pairs_Hz remains in Hz until create. Fields with conflicting native metadata are refused rather than overwriting existing inputs.

bath.ids and bath.groups retain original candidate-site labels after occupancy and exclusion; coordinates_A, isotopes, configurations, and original candidate metadata retain their finite-bath meanings. bath.spins and bath.spin_groups explicitly map retained sites/groups to physical indices in the augmented sys. Existing sys/inter own the physical model; the utility does not return simulation parameters or invoke any old public wrapper.

If neither Zeeman nor susceptibility input is present, introducing a generated Zeeman table preserves create's magnet-frequency default for pre-existing electrons using its free-electron g value. An explicitly empty Zeeman declaration remains distinct and is not replaced; scalar, eigenvalue/Euler, and matrix inputs retain their entries.

Tensor-free appended electrons also retain the standalone magnet-frequency and unit dipolar-scaling defaults when an existing Zeeman or susceptibility declaration disables create's global fallback. Only the new electron slots receive the free-electron g tensor; existing empty, scalar, eigenvalue/Euler, and matrix declarations remain unchanged, and explicitly supplied zero bath tensors remain zero.

When `inter.giant` is supplied, its native `coeff` and `euler` tables gain one empty entry per appended physical spin. Existing site entries are unchanged; `create` still rotates the supplied coefficients and converts Hz to rad/s. An absent giant-spin declaration remains absent.
