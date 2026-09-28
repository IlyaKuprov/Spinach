# kernel/utilities/corrfun.m

- Signature: `[weights,rates,states]=corrfun(spin_system,n,k,m,p,q)`

## Purpose

Computes a Wigner matrix element correlation function under isotropic, axial, or rhombic rotational diffusion.

## Physical / mathematical content

- The five indices specify the ensemble-averaged Wigner function product `<D{n}{k,m}(0)*D{n}{p,q}(t)'>`.
- Wigner function indices are sorted in descending order: for `n=2`, input indices `k=[1 2 3 4 5]` represent `[2 1 0 -1 -2]`.
- Second-rank rotational correlation times (as per Spinach input) are updated automatically if other ranks are specified.

## Numerical / algorithmic content

- The function processes each chemical species separately and identifies its basis states.
- One rotational correlation time selects the isotropic model; two select the axial model; three select the rhombic model. Diffusion coefficients are calculated as `1/(6*tau_c)` from the supplied correlation times.
- Isotropic diffusion uses one component with weight `krondelta(k,p)*krondelta(m,q)/(2*n+1)` and rate `-n*(n+1)*D`. Axial diffusion uses the same weight and a rate determined by `D_ax`, `D_eq`, `n`, and the input index `m`.
- Rhombic diffusion computes five decay rates and their weights from `Dxx`, `Dyy`, and `Dzz`. It requires `n=2` and rejects diffusion coefficients that are nearly equal.
- Input checks require the `sphten-liouv` basis formalism and real, integer Wigner indices within their permitted ranges.

## Parameters / inputs

- spin_system -the output of [[create.m]] to which ro-
- tational correlation time should have
- been supplied. For a single correlation
- time, the isotropic rotational diffusi-
- on model is used; a vector with two
- correlation times is assumed to be cor-
- relation times for rotation around and
- perpendicularly to the main axis res-
- pectively); a vector with three corre-
- lation times is assumed to be the cor-
- relation times for the rotation around
- the XX, YY and ZZ direction respecti-
- vely of the rotational diffusion tensor.
- n,k,m,p,q -the five indices found in the ensemble-
- averaged Wigner function product:
- <D{n}{k,m}(0)*D{n}{p,q}(t)'>

## Outputs

- weights -a cell array (one element for each che-
- mical species) of vectors listing the
- weights of the exponential components
- of the decays
- rates -a cell array (one element for each che-
- mical species) of vectors listing the
- decay rates (negative numbers) of the
- exponential components of the decays
- states -a cell array (one element for each che-
- mical species) of logical vectors indi-
- cating which states in the basis set
- belong to which chemical species
- Note: Wigner function indices are sorted in descending order, that is,
- k=[1 2 3 4 5] in the input represents [2 1 0 -1 -2] for n=2.
- Note: second rank rotational correlation times (as per Spinach input)
- will be updated automatically if other ranks are specified.

## Implementation structure

- Validates the basis formalism and indices, allocates per-species output cells, identifies each species’ basis states, then selects a diffusion model according to the number of correlation times.
- Source reference: <https://spindynamics.org/wiki/index.php?title=corrfun.m>