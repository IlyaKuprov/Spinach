# etc/build_isotopes.m

Source: [etc/build_isotopes.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/build_isotopes.m)

`isotopes=build_isotopes(source_file)` is an offline preparation step, not a runtime lookup fallback. It imports the audited TSV using explicit column names and types, derives nuclear g factors and missing nuclear gamma using CODATA2022 nuclear magneton and exact hbar, applies spin-zero and quadrupole selection rules, and writes an uncompressed `-v7` MAT table alongside the source. Gamma is angular frequency per field; magnetic moments are in nuclear magnetons, quadrupole moments in barns, abundances are fractions, and half-lives are seconds.

The builder checks unique labels, nuclear atomic/mass numbers, half-integer spins, abundance bounds, positive known half-lives, and row DOI attribution. It does not certify source correctness or completeness merely because those checks pass: literature adoption and state identity must be audited before changing the TSV. Direct particle gamma values are retained; particles do not receive nuclear g factors. The table contains row names, variable units/descriptions, and a schema marker. Rebuild the payload after an audited TSV edit, then validate it against the public `spin` API; ordinary numerical calls never reopen the TSV.
