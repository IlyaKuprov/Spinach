# experiments/spin_chem/rydmr_exp.m

- Signature: `answer=rydmr_exp(spin_system,parameters,H,R,K)`

## Purpose

Computes singlet-singlet RYDMR yields for a two-electron singlet over the specified magnetic fields and recombination rates, with exponential recombination (http://dx.doi.org/10.1080/00268979809483134).

## Physical / mathematical content

- The exponential recombination function is described in (http://dx.doi.org/10.1080/00268979809483134). Separates the field-dependent Zeeman term from the zero-field Hamiltonian and includes relaxation, chemical kinetics, and the source-defined exponential recombination term.
- Supports the source's Liouville- and Hilbert-space branches. The function warns that exponential recombination is built in and must not be combined with `inter.chem.rp_rates`.

## Numerical / algorithmic content

- Evaluates the yield over the requested field-rate pairs using the applicable formalism-specific calculation, with parallel loops over the pairs.
- Reshapes the calculated values to the output dimensions defined by the field and rate inputs.

## Parameters / inputs

- parameters.fields -row vector of field values, Tesla; the
- primary magnet field should be set to
- sys.magnet=1 for normalisation purposes
- parameters.rates -row vector of singlet recombination
- rate constants, Hz
- parameters.electrons -numbers identifying the two electrons
- in the isotope list, e.g. [1 2]
- parameters.needs -must contain 'zeeman_op', this is an
- instruction to the kernel to provide a
- separate Zeeman operator for field sweep
- purposes

## Outputs

- A -a matrix of singlet yields with dimensions
- matching the sizes of parameters.rates and
- parameters.fields
- Note: exponential recombination kinetics is built into this func-
- tion, do not combine with inter.chem.rp_rates parameter.

## Implementation structure

- Validates the spin formalism, operator dimensions, magnetic-field and rate inputs, and the two electron indices; the Hilbert-space branch also checks the relaxation operator.
- Constructs the two-electron singlet and Zeeman contribution, evaluates each field-rate case in the selected formalism, and returns the reshaped yield array.
