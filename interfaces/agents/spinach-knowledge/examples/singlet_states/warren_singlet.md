# examples/singlet_states/warren_singlet.m

Source: [examples/singlet_states/warren_singlet.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/singlet_states/warren_singlet.m)

## Purpose

Construct and diagonalise a full liquid-state Redfield relaxation superoperator for a two-spin 14N system. The example is intended to examine the source comment's case of long-lived states under dipolar and quadrupolar relaxation; it reports the sorted relaxation-superoperator eigenvalues rather than a propagated trajectory or an experimentally measured lifetime.

## Spin system and interactions

The field is 14.1 T. The two 14N coordinates are entered as [0.0, 0.0, 0.0] and [0.6, 0.8, 1.0]. For each nucleus, a quadrupolar coupling matrix is constructed with eeqq2nqi(1.25e6, 0.25, 1, [0, 0, 0]). The example uses Redfield relaxation, zero equilibrium, lab-frame retention, and a 5 ns correlation time. Relaxation integration and zero-value tolerances are both 1e-5.

## Calculation and output

The complete sphten-liouv basis is selected with no approximation. The script builds the spin system, obtains the relaxation matrix with relaxation(spin_system), and evaluates sort(eig(full(R))). Thus the explicit numerical output requested by the script is the ordered spectrum of the full relaxation matrix. The source comments describe the motivation as testing circumstances in which long-lived states can resist dipolar, CSA, and quadrupolar relaxation; no CSA interaction is separately specified in this script, and the source-only analysis does not establish a lifetime or numerical result.

## Citation

The source file names Warren and Kuprov in its comments but supplies no DOI or publication citation.
