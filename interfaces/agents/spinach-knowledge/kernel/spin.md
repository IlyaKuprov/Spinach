# kernel/spin.m

- Signature: `[gamma,multiplicity,data]=spin(name)`
- Source: [kernel/spin.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/spin.m)
- Wiki: [spin.m](https://spindynamics.org/wiki/index.php?title=spin.m)

## Numerical and metadata contracts

A character row vector selects a canonical physical label. Existing one-output calls return signed angular magnetogyric ratio in rad/(s*T); two-output calls additionally return `2*I+1`. Divide gamma by `2*pi` for Hz/T. The third output is a table row containing spin, abundance, quadrupole moment, radioactive half-life in seconds, source qualifications, and field-labelled academic DOIs. `[~,~,data]=spin('table')` returns the full physical table; the `table` selector requires three outputs. Returned tables have value semantics: editing them does not modify the persistent database.

Unknown labels raise `spin:unknown_isotope`. A known row without a confirmed spin or signed gamma raises `spin:data_unavailable` for numerical calls. Such rows remain inspectable through the full table. Tentative/systematic spins, unsigned experimental moments, unknown abundance, and uncertain half-lives are not silently made into confirmed simulation inputs. Metadata does not automatically introduce abundance weighting or radioactive decay into a simulation.

## Labels and special cases

Nuclides use mass number plus element symbol, such as `13C` or `195Pt`; separate nuclear states have explicit suffixes, such as `99Tc_m`. `180Ta` is the short-lived spin-one ground state; natural tantalum-180 abundance is assigned to `180Ta_m`, the spin-nine isomer. Unknown state energies and identities remain qualified rather than being forced into a ground-state assignment.

Particle rows include `E`/`E+`, `N`/`antiN`, `M`/`M+`, `anti1H`, and static-moment hyperons and their explicit antiparticles. `1H` is the proton's existing nuclide row, not a duplicate particle. Antiparticle properties inferred under CPT are labelled as assumptions, not independent measurements. Magnetic moments use nuclear magnetons for both row kinds; nuclear g factors are inapplicable to particles.

`G` is a ghost: gamma zero, multiplicity one. `E#` denotes a high-spin electron with integer multiplicity at least two; `C#`, `V#`, and `T#` denote cavity, phonon, and transmon modes with integer level count at least three and gamma zero. These abstract specifications have typed empty metadata, not physical table rows.

## Units, sources, and retrieval

The authoritative editable dataset is `etc/isotopes.tsv`; `etc/build_isotopes.m` creates its typed, uncompressed MAT release projection. Source selection and qualifications are documented in `etc/isotopes_sources.md`. Spin-zero magnetic moments and gamma are zero; quadrupole moments vanish for spin below one. Missing properties are NaN, not invented zeros. Interval-only natural abundance has no scalar midpoint. Stable half-life Inf means classified stable/no observed decay, not a measured infinite lifetime; alternative environmental lifetimes remain explicit in notes.

Nuclear states whose adopted half-life is below one second (value, estimate, or upper limit) are excluded. Stable nuclei, exactly one second, unknown lifetimes, and unresolved lower limits remain; particle rows are not subject to the isotope cutoff. Alternative environmental lifetimes in notes do not override the adopted value. Excluded labels raise the ordinary unknown-isotope error.

The MAT table is loaded once per process. Warm numeric calls use a label dictionary and cached numeric projections, without filesystem access or table slicing. Metadata rows are constructed only when requested. No online dependency, TSV fallback, or automatic regeneration is used; `clear spin` or a new process loads an updated release payload.
