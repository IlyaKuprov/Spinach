# etc/estimators/guess_csa_pro.m

- Signature: `CSAs=guess_csa_pro(aa_nums,pdb_ids,coords,options)`

## Purpose

Guesses a reasonable amide bond 15N CSA tensor anisotropy, and a reasonable 13C=O tensor anisotropy, given a local protein geome- try. The tensors are oriented roughly according to Proline is not currently handled, amino acids are assumed to be numbered from N-terminus to C-terminus. Syntax: CSAs=guess_csa_pro(aa_nums,pdb_ids,coords,options)

## Physical / mathematical content

- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

## Parameters / inputs

- aa_nums -a vector of amino acid numbers
- pdb_ids -a cell array of PDB atom identifiers
- coords -a cell array of coordinate vectors
- options.nh_csa -'tcb' (default) for Tjandra, Curtis,
- and Bodenhausen, 'bax' for Cornilescu
- and Bax, and 'pol' for Case, Polenova
- and Gronenborn eigenvalues and orien-
- tations of the CSA tensors

## Outputs

- CSAa -a cell array of 3x3 CSA tensors in ppm
- Note: these CSAs are very approximate. For accurate relaxation
- analysis you must supply your own tensors.
- Note: this is an auxiliary function that is called by protein.m
- protein import module. Direct calls are discouraged.

## Implementation structure

- Guesses a reasonable amide bond 15N CSA tensor anisotropy, and a
- reasonable 13C=O tensor anisotropy, given a local protein geome-
- try. The tensors are oriented roughly according to
- Proline is not currently handled, amino acids are assumed to be
- numbered from N-terminus to C-terminus. Syntax:
- CSAs=guess_csa_pro(aa_nums,pdb_ids,coords,options)
- aa_nums -a vector of amino acid numbers
- pdb_ids -a cell array of PDB atom identifiers
- coords -a cell array of coordinate vectors
- options.nh_csa -'tcb' (default) for Tjandra, Curtis,
- and Bodenhausen, 'bax' for Cornilescu
- and Bax, and 'pol' for Case, Polenova
