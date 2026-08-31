# etc/estimators/guess_csa_pro.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/estimators/guess_csa_pro.m`
- Signature: `CSAs=guess_csa_pro(aa_nums,pdb_ids,coords,options)`
- Total lines: 338

## Purpose

Guesses a reasonable amide bond 15N CSA tensor anisotropy, and a reasonable 13C=O tensor anisotropy, given a local protein geome- try. The tensors are oriented roughly according to Proline is not currently handled, amino acids are assumed to be numbered from N-terminus to C-terminus. Syntax: CSAs=guess_csa_pro(aa_nums,pdb_ids,coords,options)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 44-45: Set default peptide bond CSA options; implemented by `if ~isfield(options,'nh_csa'), options.nh_csa='tcb'; end`.
- Lines 47-48: Check consistency; implemented by `grumble(aa_nums,pdb_ids,coords,options)`.
- Lines 50-51: Preallocate CSA array; implemented by `CSAs=cell(numel(pdb_ids),1)`.
- Lines 53-54: Number the atoms; implemented by `numbers=1:numel(coords)`.
- Lines 56-57: Loop over amino acids; implemented by `for n=1:(max(aa_nums)-1)`.
- Lines 59-62: Assign amide nitrogen CSAs; implemented by `if ismember('C',pdb_ids(aa_nums==n))&& ismember('N',pdb_ids(aa_nums==(n+1)))&& ismember('H',pdb_ids(aa_nums==(n+1)))`.
- Lines 64-65: Get C coordinates; implemented by `local_coords=coords(aa_nums==n)`.
- Lines 68-69: Get N coordinates; implemented by `local_coords=coords(aa_nums==(n+1))`.
- Lines 72-73: Get H coordinates; implemented by `local_coords=coords(aa_nums==(n+1))`.
- Lines 76-77: Get the primary directions; implemented by `N_CO_vec=C-N; N_H_vec=H-N`.
- Lines 79-80: Double-check the distances; implemented by `if (norm(N_CO_vec,2)>2.0)||(norm(N_H_vec,2)>2.0)`.
- Lines 84-85: Make ZZ eigenvector collinear with N-CO bond; implemented by `zz_eigvec=N_CO_vec`.
- Lines 88-89: Make YY eigenvector perpendicular to the peptide plane; implemented by `yy_eigvec=cross(N_CO_vec,N_H_vec)`.
- Lines 92-93: Make XX eigenvector perpendicular to the other two; implemented by `xx_eigvec=cross(yy_eigvec,zz_eigvec)`.
- Lines 96-97: Build the eigenvalue matrix; implemented by `switch options.nh_csa`.
- Lines 99-100: Cornilescu and Bax; implemented by `case 'bax'`.
- Lines 105-106: Tjandra, Curtis, and Bodenhausen; implemented by `case 'tcb'`.
- Lines 113-114: Case, Polenova, Gronenborn; implemented by `case 'pol'`.

### Control flow inferred from the code

- Line 45: conditional branch on `~isfield(options,'nh_csa'), options.nh_csa='tcb'; end`.
- Line 57: `for` loop over `n=1:(max(aa_nums)-1)`.
- Line 60: conditional branch on `ismember('C',pdb_ids(aa_nums==n))&&`.
- Line 80: conditional branch on `(norm(N_CO_vec,2)>2.0)||(norm(N_H_vec,2)>2.0)`.
- Line 97: dispatches on `options.nh_csa`; cases `'bax'`, `'tcb'`, `'pol'`.
- Line 149: `for` loop over `n=1:(max(aa_nums)-1)`.
- Line 152: conditional branch on `ismember('C',pdb_ids(aa_nums==n))&&`.
- Line 172: conditional branch on `(norm(N_C_vec,2)>2.0)||(norm(C_CA_vec,2)>2.0)`.
- Line 214: `for` loop over `n=2:(max(aa_nums)-1)`.
- Line 217: conditional branch on `ismember('H',pdb_ids(aa_nums==n))&&`.
- Line 242: conditional branch on `(norm(N_H_vec,2)>1.2)||(norm(N_CA_vec,2)>1.6)||(norm(N_C_vec,2)>1.5)`.
- Line 259: dispatches on `options.nh_csa`; cases `'bax'`, `'tcb'`, `'pol'`.

### Key state/data transformations

- Lines 51: computes `CSAs` using `CSAs=cell(numel(pdb_ids),1)`.
- Lines 54: computes `numbers` using `numbers=1:numel(coords)`.
- Lines 77: computes `N_CO_vec` using `N_CO_vec=C-N; N_H_vec=H-N`.
- Lines 85: computes `zz_eigvec` using `zz_eigvec=N_CO_vec`.
- Lines 89: computes `yy_eigvec` using `yy_eigvec=cross(N_CO_vec,N_H_vec)`.
- Lines 93: computes `xx_eigvec` using `xx_eigvec=cross(yy_eigvec,zz_eigvec)`.
- Lines 103: computes `D` using `D=diag([-108.0 62.0 46.0])`.
- Lines 127: computes `V` using `V=[xx_eigvec yy_eigvec zz_eigvec]`.
- Lines 134: computes `CSAs{nitrogen_number}` using `CSAs{nitrogen_number}=V*D*V'`.
- Lines 169: computes `N_C_vec` using `N_C_vec=C-N; C_CA_vec=CA-C`.
- Lines 199: computes `CSAs{carbon_number}` using `CSAs{carbon_number}=V*D*V'`.
- Lines 239: computes `N_H_vec` using `N_H_vec=N-H; N_CA_vec=N-CA; N_C_vec=N-C; H_CA_vec=H-CA`.
- Lines 296: computes `CSAs{proton_number}` using `CSAs{proton_number}=V*D*V'`.

### Local helper functions

- Line 313: `grumble()` — `function grumble(aa_nums,pdb_ids,coords,options)`.
  - Representative operation: `if ~isnumeric(aa_nums)`.
  - Representative operation: `error('aa_nums must be a vector of integers.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isfield()`, `grumble()`, `ismember()`, `pdb_ids()`, `coords()`, `strcmp()`, `cross()`, `numbers()`, `local_numbers()`, `num2str()`, `iscell()`, `isstruct()`, `ischar()`.
