# kernel/states/deut_pair.m

- Signature: [S,T,Q,Tc,Qc] = deut_pair(spin_system,spin_a,spin_b,options)
- Source: [kernel/states/deut_pair.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/states/deut_pair.m)
- Wiki: [deut_pair.m](https://spindynamics.org/wiki/index.php?title=deut_pair.m)
- Basis construction reference: [Chemistry Physics Letters 296, 179 (1998)](https://doi.org/10.1016/S0009-2614(98)00784-2)

## Purpose

Builds population and coherence state combinations for a pair of selected spin-1 sites. The implementation requires multiplicity 3 for each selected spin and supports the 'zeeman-hilb', 'zeeman-liouv', and 'sphten-liouv' formalisms.

## Basis and formalism

The single-spin coordinate basis is ordered as alpha = [1;0;0], beta = [0;1;0], gamma = [0;0;1]. Every pair product is kron(first,second): spin_a is the first factor and spin_b the second. The helper uses the tensor products of irr_sph_ten(3) elements as its product basis and normalises each product basis operator to Frobenius norm 1. Spinach state descriptors are formed as T followed by the L,M indices for each selected site.

The source constructs these labelled population states:

- S0 = (alpha x gamma - beta x beta + gamma x alpha)/sqrt(3).
- T order is [T+,T0,T-]: T+ = (alpha x beta - beta x alpha)/sqrt(2); T0 = (alpha x gamma - gamma x alpha)/sqrt(2); T- = (beta x gamma - gamma x beta)/sqrt(2).
- Q order is [Q++,Q+,Q0,Q-,Q--]: Q++ = alpha x alpha; Q+ = (alpha x beta + beta x alpha)/sqrt(2); Q0 = (alpha x gamma + 2 beta x beta + gamma x alpha)/sqrt(6); Q- = (beta x gamma + gamma x beta)/sqrt(2); Q-- = gamma x gamma.

For 'zeeman-hilb', each generated state rho is normalised by its Frobenius norm. For 'zeeman-liouv' and 'sphten-liouv', the full state is normalised by its 2-norm. Other formalisms raise an unsupported-formalism error. The source cautions that the labelled states are not irreducible spherical tensors.

## Inputs

- spin_system — compiled Spinach system, including the formalism and selected-site multiplicities.
- spin_a, spin_b — indices of the two selected sites; both must have multiplicity 3.
- options.dephasing — numeric scalar 0 or 1. If the options argument is omitted, the helper sets it to 0. With 1, it skips terms unless both component indices M1 and M2 are zero; with 0, it retains all terms. Supplying options without this field is rejected.

## Outputs

- S — the singlet-labelled population combination.
- T — three triplet-labelled population combinations in the order [T+,T0,T-].
- Q — five Q-labelled population combinations in the order [Q++,Q+,Q0,Q-,Q--].
- Tc — four triplet coherence combinations, in source order: T0 to T-, T+ to T0, T- to T0, T0 to T+.
- Qc — eight Q-state coherence combinations, in source order: Q- to Q--, Q0 to Q-, Q+ to Q0, Q++ to Q+, Q-- to Q-, Q- to Q0, Q0 to Q+, Q+ to Q++.

The returned population and coherence entries are Spinach states in the selected formalism, not scalar probabilities. Coherence lists use the labels and ordering in the source.
