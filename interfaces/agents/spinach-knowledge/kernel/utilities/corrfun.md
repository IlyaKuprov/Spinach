# kernel/utilities/corrfun.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/corrfun.m`
- Signature: `[weights,rates,states]=corrfun(spin_system,n,k,m,p,q)`
- Total lines: 162

## Purpose

Wigner matrix element correlation function under isotropic, axial, and rhombic rotational diffusion. Syntax: [weights,rates,states]=corrfun(spin_system,n,k,m,p,q)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

- Wigner matrix element correlation function under isotropic, axial,
- and rhombic rotational diffusion. Syntax:
- [weights,rates,states]=corrfun(spin_system,n,k,m,p,q)
- spin_system -the output of [[create.m]] to which ro-
- tational correlation time should have
- been supplied. For a single correlation
- time, the isotropic rotational diffusi-
- on model is used; a vector with two
- correlation times is assumed to be cor-
- relation times for rotation around and
- perpendicularly to the main axis res-
- pectively); a vector with three corre-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `krondelta()`, `strcmp()`, `any()`.
