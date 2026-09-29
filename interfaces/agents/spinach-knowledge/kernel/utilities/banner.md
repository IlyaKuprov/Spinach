# kernel/utilities/banner.m

## Purpose

Prints console banners for the Spinach kernel. This is an internal kernel function; user calls are discouraged.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/banner.m>

## Behaviour

- Syntax: `banner(spin_system,identifier)`.
- The function first validates `identifier` via the internal `grumble` helper, which errors with `'identifier must be a character string.'` if `identifier` is not a character array.
- A `switch` on `identifier` selects one of six banner cases, each printing a bordered banner block through `report(spin_system,...)` calls:
  - `'version_banner'` — prints the SPINACH v2.13 banner, including hyperlinked lines for the author list (`https://spindynamics.org/wiki/index.php?title=Spinach_developer_team`), documentation (`https://spindynamics.org/wiki/index.php?title=Main_Page`), and book (`https://link.springer.com/book/10.1007/978-3-031-05607-9`), plus an `MIT License` line.
  - `'spin_system_banner'` — prints a `SPIN SYSTEM` banner.
  - `'basis_banner'` — prints a `BASIS SET` banner.
  - `'sequence_banner'` — prints a `PULSE SEQUENCE` banner.
  - `'optimcon'` — prints an `OPTIMAL CONTROL` banner.
  - `'optimisation'` — prints an `OPTIMISATION` banner.
- Any other identifier triggers `error('unknown banner.')`.
- All banner blocks are delimited by lines of `=` characters with blank `report(spin_system,' ')` lines before and after.

## Inputs and outputs

- `spin_system` — spin system object passed through to `report` for console output.
- `identifier` — character string selecting the banner; one of `'version_banner'`, `'spin_system_banner'`, `'basis_banner'`, `'sequence_banner'`, `'optimcon'`, `'optimisation'`.
- Returns nothing; output is console text via `report`.

## References

- Spinach Wiki page for this function: <https://spindynamics.org/wiki/index.php?title=banner.m>
- Author list: <https://spindynamics.org/wiki/index.php?title=Spinach_developer_team>
- Documentation: <https://spindynamics.org/wiki/index.php?title=Main_Page>
- Book: <https://link.springer.com/book/10.1007/978-3-031-05607-9>
