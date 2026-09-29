# kernel/spin.m

- Signature: [gamma,multiplicity] = spin(name)
- Source: [kernel/spin.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/spin.m)
- Wiki: [spin.m](https://spindynamics.org/wiki/index.php?title=spin.m)

## Purpose

Looks up the magnetogyric ratio and level multiplicity for a named isotope or supported particle. The input name identifies the database entry; the outputs are not a full isotope-metadata record.

## Meaning of outputs

- gamma is the magnetogyric ratio in rad/(s*T), an angular-frequency-per-field unit. It is not in Hz/T; divide by 2*pi to express the same value in Hz/T.
- multiplicity is the number of energy or population levels stored or assigned by this function. It is not the nuclear spin quantum number I, and spin does not return I as an output.
- The function has exactly two outputs: gamma first and multiplicity second. Requesting one output returns gamma.

## Lookup cases

The input can be an isotope name such as '15N' or '195Pt'. The source also handles 'G' (ghost spin: gamma = 0, multiplicity = 1), 'N' (neutron), and 'M' (muon). Parameterised names are 'E#' for a high-spin electron, 'C#' for an electromagnetic cavity mode, 'V#' for a phonon mode, and 'T#' for a transmon; # supplies the multiplicity or level count. The electron case requires at least two levels. Cavity, phonon, and transmon cases require at least three levels and assign gamma = 0.

The isotope cases assign gamma and multiplicity from the table in the source. An unrecognised name raises an unknown-isotope error. The function returns these values; it does not print a summary.

## Notes

The source marks some values as needing verification and says entries without a stated source were taken from Google and should be double-checked before production calculations. These are database comments, not independent verification of isotope data.
