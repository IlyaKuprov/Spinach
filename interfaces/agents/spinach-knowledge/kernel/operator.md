# kernel/operator.m

- Signature: `A=operator(spin_system,operators,spins,operator_type,format)`

## Meaning

Builds the requested spin operator or Liouville-space superoperator in the basis selected by `spin_system.bas.formalism`. The accepted specification forms are paired character strings (operator labels plus isotope/type), a character label plus a numeric row of spin indices, or paired cell arrays of labels and numeric spin indices. The first two forms sum the corresponding single-spin operators; paired cells specify a product. Source examples are `operators='Lz'; spins='13C'`, `operators='Lz'; spins=[1 2 4]`, and `operators={'Lz','L+'}; spins={1,2}`.

Documented operator labels include `E`, `Lz`, `Lx`, `Ly`, `L+`, `L-`, `Tl,m`, and Zeeman-basis central-transition labels `CTx`, `CTy`, `CTz`, `CT+`, `CT-`. String spin selectors may be isotope names or `electrons`, `nuclei`, or `all`.

## Basis and index construction

The routine passes the request to `human2opspec`, then assembles each returned specification. In Zeeman formalisms it maps each per-spin linear specification through `lin2lm` to `L,M`, selects `irr_sph_ten(mults(k),L){L-M+1}` for that spin, and Kronecker-products the local matrices in specification order. In `zeeman-liouv` the resulting Hilbert operator is converted by `hilb2liouv`. In `sphten-liouv` it calls `superop` for each specification and confines any identity term to the selected substance. The code uses a sum of the resulting terms with their returned coefficients. The substance selector is defined before parallel assembly even for single-substance sum requests, where no block confinement is needed.

For Liouville calculations the documented `operator_type` values are `left`, `right`, `comm` (default), and `acomm`. Hilbert formalisms ignore this option. A Liouville product request represents the superoperator of the full product operator, not a product of separately generated single-spin superoperators.

## Output and guards

- `format='csc'` (default) returns a sparse square matrix. Its dimension is `bas.offsets(end)`, the compiled direct-sum dimension.
- `format='xyz'` returns the nonzero entries as `[row,column,value]` triplets, with row and column indices from MATLAB `find`.
- It requires basis information, nonempty spin selections, and unique positive integer spin indices; paired cell arrays must have matching lengths and character operator labels. The format must be `csc` or `xyz`. The local guard requires `operator_type` to be a character string; the downstream superoperator/conversion routines handle its meaning.
- Product specifications must lie within one substance, including explicitly specified identity factors. Isotope and numeric-vector requests remain sums of single-spin operators. Zeeman tensor products are built only from the hosting substance's spins and offset into its Hilbert or Liouville block; other blocks remain zero. Isotope and numeric selections sum these local operators without inter-species coherences.
- Explicit `E` and `T0,0` requests act only in the selected block: left/right actions give its identity, anticommutation gives twice that identity, and commutation gives zero. Local products contribute once; numeric and isotope selections retain one contribution per selected spin. This applies to both output formats.
- Caching is optional: `op_cache` in `spin_system.sys.enable` enables the cache path when its ValueStore is available.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operator.m)
- [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=operator.m)
