# interfaces/pdb_bmrb/protein.m

- Signature: `[sys,inter,aux]=protein(pdb_file,bmrb_file,options)`

## Purpose

Import protein data from PDB and BMRB files into Spinach structures. The function matches atom coordinates and chemical-shift assignments, estimates J-couplings and chemical-shift anisotropies (CSAs), and returns system, interaction, and auxiliary data.

## Physical / mathematical content

- PDB data provide residue and atom identifiers and coordinates; BMRB data provide chemical-shift assignments.
- The function estimates scalar J-couplings and stores CSA tensors as 3x3 matrices. The output records chemical shifts in ppm and scalar couplings in Hz.

## Numerical / algorithmic content

- Reads the selected PDB molecule and BMRB file, then matches assignments by residue number and atom identifier. A residue-type mismatch between the files is an error.
- Uses atom-name heuristics to replicate some BMRB shifts for corresponding protons and symmetric aromatic-ring positions. Some unassigned terminal or side-chain groups are ignored.
- Calls `guess_j_pro` and `guess_csa_pro` to estimate scalar couplings and CSAs, then applies atom selection and the missing-shift policy. Isotopes are assigned as `1H`, `13C`, or `15N` from atom identifiers.
- When deuteration is requested, specified protons are converted to `2H` and their labels use a `D` prefix. Their shifts are retained, and their scalar couplings are scaled by `spin('2H')/spin('1H')`.

## Parameters / inputs

- `pdb_file`: string naming the PDB file.
- `bmrb_file`: string naming the BMRB file.
- `options.select`: `'backbone'` imports the backbone through CB and HB; `'backbone-minimal'` imports only the backbone; `'backbone-hsqc'` also includes GLN and ASN side-chain amide groups; `'all'` imports every atom assigned in BMRB. Alternatively, a list of PDB atom serial numbers can be supplied; every number must be present in the file. Unsupported types (oxygen, sulphur, and OH protons) are dropped with a warning. Atoms without a BMRB assignment are kept or deleted according to `options.noshift`.
- `options.pdb_mol`: molecule number to read when the PDB file contains multiple molecules.
- `options.noshift`: `'keep'` places unassigned atoms between -1 and 0 ppm; `'delete'` removes them from the system.
- `options.deuterate`: a cell array of character strings naming PDB atom identifiers for protons to replace with deuterons, or `'non-Me'` to deuterate everything except methyl groups.
- `options.nh_csa`: peptide-bond CSA choice. `'bax'` gives H [6.00 0.00 -6.00] and N [-108.0 62.0 46.0] ppm; `'tcb'` gives H [7.00 0.00 -7.00] and N [-125.0 45.0 80.0] ppm; `'pol'` gives H [6.66 0.66 -7.33] and N [-92.4 34.7 57.7] ppm. The documented default is `'tcb'`.

If `options` is omitted, the function sets `select='all'`, `pdb_mol=1`, `noshift='keep'`, and `deuterate={}`. If supplied, the options structure must include `select`, `pdb_mol`, and `noshift`; an omitted `deuterate` field defaults to `{}`.

## Outputs

- `sys.isotopes`: `Nspins x 1` cell array of isotope strings.
- `sys.labels`: `Nspins x 1` cell array of standard IUPAC protein atom labels.
- `inter.coordinates`: `Nspins x 3` matrix of coordinates in angstroms.
- `inter.zeeman.scalar`: `Nspins x 1` cell array of isotropic chemical shifts in ppm. The implementation assigns this field; the source header instead names `inter.zeeman.iso`.
- `inter.zeeman.matrix`: `Nspins x 1` cell array of `3x3` CSA matrices in ppm.
- `inter.coupling.scalar`: `Nspins x Nspins` cell array of scalar couplings in Hz.
- `aux.pdb_aa_num`: PDB amino-acid number for each spin.
- `aux.pdb_aa_typ`: PDB amino-acid type for each spin.
