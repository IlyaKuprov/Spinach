# examples/spin_chemistry/singlet_yield_nz_1.m

Source: [examples/spin_chemistry/singlet_yield_nz_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/spin_chemistry/singlet_yield_nz_1.m)

## Purpose and cited context

Compare Redfield relaxation with the lifetime-shifted Nakajima-Zwanzig (NZ) kernel for singlet recombination yield in a triplet-born benzophenone-ketyl/thiyl radical pair. Source comments describe a viscous ionic liquid and say the parameters are representative of TMPA-TFSA measurements by the Wakasa group. The two DOIs are citations in the source comments, not evidence that these calculations reproduce or independently verify experimental results: [10.1021/jp074331a](https://doi.org/10.1021/jp074331a) and [10.1039/c2cp23747d](https://doi.org/10.1039/c2cp23747d).

The source comment characterises the cage lifetime as tens of nanoseconds and the rotational correlation time as a few nanoseconds, with k*tau_c near 0.1. The code uses k_cage=3e7 Hz and tau_c=3.3e-9 s, whose product is 0.099.

## Spin model and relaxation

For each calculation the helper sets a two-electron/two-proton system (E, E, 1H, 1H) and electron scalar Zeeman values 2.0032 and 2.0080. The ketyl-proton principal hyperfine inputs are mt2hz([0.14 0.14 0.59]) and mt2hz([0.21 0.21 0.66]); the latter is rotated by Euler angles [pi/7 pi/5 pi/3]. These numbers are supplied to mt2hz as written; the source does not state a separate unit label for them.

Cage recombination is one explicit unselected first-order loss record at total rate k_cage; the singlet-yield prefactor remains k_cage/2. The relaxation setting is selected by theory, with equilibrium zero, lab-frame retention, rlx_dfs='keep', and tau_c={tau_c}. For the NZ branch, the source additionally selects nz_shift=k_cage and nz_onshell=false. The basis is the full sphten-liouv basis (bas.approximation='none').

The initial operator is constructed explicitly as rho0=(EE-S)/3, where S is the singlet projector for spins 1 and 2 and EE is the two-electron identity state. The script builds Hamiltonian, relaxation, and kinetics superoperators and solves the reaction integral with bicg on L=H+1i*R+1i*K. It then evaluates a normalised fractional singlet yield from the singlet and triplet operators and the cage rate. Thus the initial spin preparation and singlet-yield observable are explicit in this example rather than inferred from a plotting routine.

## Parameter sweeps and plots

The main field grid contains 15 points from 0.5e-3 to 50e-3 T, displayed in mT. At each field the code evaluates both theories and prints a table headed with field in T and the Redfield and NZ yields. A second sweep uses drain rates [1e6 3e6 1e7 3e7 1e8] Hz at 3e-3 T, reporting both yields and plotting NZ minus Redfield yield against tau_c*drains. The figures visualise calculated outputs; the source file itself does not provide a measured yield dataset.

## Linear-solver accuracy

The shipped `bicg` relative residual target is 1e-8. Comparing the full five-rate ladder before and after reaction-record migration gives relative yield differences of 4.52030e-10 (Redfield) and 4.19909e-10 (NZ). Changing only that target to 1e-12 in both calculations reduces these differences to 1.03858e-14 and 5.50697e-15, respectively. This measured tolerance refinement identifies the ladder residual as linear-solver accuracy rather than a resolved physical change; the example retains its original solver target. The complete 15-field curves at the shipped target differ by 1.39748e-11 and 4.50731e-11.
