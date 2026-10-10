# kernel/hamiltonian.m

- Signature: `[I,Q]=hamiltonian(spin_system,operator_type)`
- Default: `operator_type='comm'`

## Purpose and outputs

Builds the Hamiltonian terms represented by a configured Spinach spin system. `I` contains the rotationally invariant terms; when a second output is requested, `Q` contains irreducible components of the anisotropic spin Hamiltonian. It is a cell array by rank; `Q{r}` is a `(2r+1)-by-(2r+1)` cell array indexed by operator and coefficient projections, and `Q{r}{k,m}` stores the corresponding basis-space matrix. Use `orientation.m` to combine those components at a chosen spin-system orientation. A one-output call skips construction of `Q` and the anisotropic blocks guarded by `nargout>1`.

In the sphten direct sum, interaction operators are placed in their hosting substance blocks. The parallel construction retains `chem.parts` so global interaction indices can be translated into local descriptor columns.

The routine uses the active basis and formalism in `spin_system`; it does not select a basis. It requires basis information and the interaction-strength assumptions established by `assume()`.

## Interactions handled

- Zeeman terms read `inter.zeeman.matrix` and `inter.zeeman.strength` for each spin. Isotropic handling accepts `full`, `secular`, and `ignore`. The secular value subtracts that spin's base frequency; full retains it. For anisotropy, `full` builds rank-1 and rank-2 components, `secular` retains the rank-2 zero-projection `Lz` term, and `ignore` omits the anisotropic term.
- When `Q` is requested, significant same-spin diagonal entries of `inter.coupling.matrix` are processed as nuclear-quadrupole-interaction terms, with `strong`, `secular`, or `ignore` strength assumptions. The strong form builds all five rank-2 projections; the secular form retains the zero-projection term. Quadrupolar nuclei have spin greater than 1/2; this function uses the interaction setup supplied in `spin_system` rather than applying a spin-number guard here.
- Distinct-spin entries of `inter.coupling.matrix` supply isotropic and anisotropic bilinear pair couplings. The dispatch accepts `strong`, `secular`, `z*`, `*z`, `zz`, `z+`, `z-`, `+z`, `-z`, `T(L,+1)`, `T(L,-1)`, and `ignore`. For the isotropic part, `strong` and `secular` use `Lz*Sz` plus half-weighted `L+*S-` and `L-*S+`; `z*`, `*z`, and `zz` retain `Lz*Sz`, while the anisotropic-only labels and `ignore` omit the isotropic part. The function treats the supplied matrices as generic coupling tensors; it does not infer a chemical or physical label from them.
- `inter.giant.coeff` supplies per-spin spherical-tensor terms of the declared ranks; `inter.giant.strength` selects `strong` (all projections), `secular` (zero projection), or `ignore` handling. These anisotropic terms are built only when `Q` is requested.
- If `inter.modes` is present, the mode block requires `inter.modes.strength` (otherwise it errors and directs the caller to `assume()`). It handles frequency assumptions `full`, `offset`, or `ignore`; anharmonicity `full` or `ignore`; exchange `strong`, `rwa`, `nonelec`, or `ignore`; and Kerr, longitudinal, dispersive, `coupling_mod`, and `zeeman_mod` terms with `full` or `ignore` assumptions. The `strong` exchange branch retains counter-rotating terms, `rwa` omits them, and `nonelec` skips electron–mode pairs. With `offset` frequencies, connected mode/spin couplings determine carrier references; conflicting spin or mode carriers error, while an unreferenced nonzero-frequency mode is reported at its laboratory frequency. Only modes whose `comp.types` entry is `C`, `V`, or `T` are included. As stated in the source header, bosonic-mode terms enter `I` at their input orientation; `Q` and orientation handling apply to the spin subsystem, not mode rotations.

## Formalism and operator effects

`operator_type` is checked against `left`, `right`, `comm`, and `acomm`. In Hilbert space the type is ignored by operator construction. In Liouville space, using column-wise vectorisation of a state matrix `X` and Hamiltonian term `H`, the generated superoperators act as follows:

- `left`: `H X`
- `right`: `X H`
- `comm`: `H X - X H`
- `acomm`: `H X + X H`

Thus `comm` uses `[H,X]`, not `[X,H]`. The function constructs these operators; it does not apply a time-evolution sign or propagate a state.

## Coefficients, thresholds, and term indexing

The coefficients are assembled with the explicit signs and factors in the interaction descriptors. The source displays reported coefficients in `Hz` by dividing by `2π`; the operator assembly uses the stored angular-frequency coefficients and does not add a global `2π` or reverse their signs. Pair tensors are first screened by norm against `2π * inter_cutoff`. For descriptor rows, the code removes a row only when its isotropic coefficient is below `liouv_zero` and either the absolute-sum of `irr_comp` or of `ist_coeff` is also below `liouv_zero`.

The function packs eight descriptor slots per spin into `D1` and a 3-by-3 operator grid per selected spin pair into `D2`, filters insignificant rows, concatenates `descr=[D1;D2]`, and sets `nterms=size(descr,1)`. It shuffles the descriptor rows with `randperm` before construction. The `parfor` loop runs `n=1:nterms`; iteration `n` reads `descr(n,:)` and stores that row's one-spin or two-spin operator in `oper{n}`. The descriptor's `nL/nS` and `opL/opS` fields map to spin indices and local operator specifications; `isotropic`, `ist_coeff`, and `irr_comp` map its scalar and irreducible coefficients. The later assembly uses those same descriptor rows to accumulate `I` and the indexed `Q` components. `I` is assembled as a sparse square matrix in the active basis; each `Q{r}{k,m}` entry is a square matrix in that basis.

When `ham_cache` is enabled, the cache key includes the descriptor, `operator_type`, isotope hash, basis hash, and giant-spin terms; if modes exist it also includes the mode data and base frequencies. A client without an existing parallel pool skips the cache; workers use their current `ValueStore`. During parallel construction the client installs a `DataQueue` progress callback, reports completed operator counts and throughput after more than five seconds or at completion, and records construction/assembly timings and average worker transfer sizes; workers do not run the client diagnostic callback.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/hamiltonian.m)
- [Spinach Wiki: hamiltonian.m](https://spindynamics.org/wiki/index.php?title=hamiltonian.m)
