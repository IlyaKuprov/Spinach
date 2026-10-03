# interfaces/jcamp/jcamp_export.m

## Purpose and syntax

`text=jcamp_export(data)` exports NMR and EPR/EMR data using a single documented
structure. The character-row result contains JCAMP-DX 5.01; `data.filename`,
when present, also receives the complete file. The [input contract and examples](../../../../jcamp/README.md)
cover file ownership, block metadata, traces, general pages, and peak tables.

## Data and numerical behaviour

Real traces select incremental XYDATA only for an exactly regular axis with
a finite initial ordinate; otherwise XYPOINTS preserves the explicit coordinates. Complex traces become
NTUPLES with distinct real and imaginary pages. General variable/page inputs
represent multidimensional coordinates, ragged sampling, and additional
quadrature or receiver channels without guessing dimension order. Explicit
pairs use the NTUPLES PROFILE display method. Shared attributes cannot force a different sampling grid onto an individual page.
Repeated coordinate sets retain distinct tables in input order; an explicit
PAGE coordinate can label replicate identity for coordinate-keyed readers.
Peaks can contain heights, widths, NMR multiplicities, and assignments. Multiple
independent datasets are enclosed in a LINK block.

Floating-point column inputs retain their order and amplitudes. AFFN uses
17 significant digits and unit scale factors; there is no FFT, rescaling,
normalisation, or unit conversion. Missing observations use `?`. The result
uses ASCII, CRLF, and lines no longer than 80 characters; an existing file is
replaced only after a complete successful write.

## Metadata and limits

NMR observation frequencies are in MHz; EMR microwave frequencies are in Hz.
EPR uses `EMR SIMULATION` or `EMR MEASUREMENT`, with explicit detection mode
and a core method identifier. Required identifiers, NMR delay pairs, and
tabulated axis units are checked, but the caller supplies the full acquisition, sample, reference, and simulation metadata relevant to the
experiment. Generated structural records cannot be overridden by aliases.
NMR ordinate units include MAGNITUDE and POWER under JCAMP 5.01; supplying
a unit label does not calculate the corresponding transform. Embedded label
periods are supported in the private namespace, not ordinary reserved labels.
The exporter does not provide integer ASDF compression, molecular structure
export, or universal vendor compatibility; JCAMP readers differ in their
support for multidimensional and irregular NTUPLES.

EMR simulation source and parameter descriptions are non-empty ASCII text,
including multiline cell vectors; numeric arrays are not descriptions.

Complex NMR quadratures use `ARBITRARY UNITS`; magnitude/power data must
be real and already transformed. NMR assignments may omit heights and method
comments; widths require a method description.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/jcamp/jcamp_export.m)
- [NMR protocol](https://iupac.org/wp-content/uploads/2021/08/JCAMP-DX_NMR_1993.pdf)
- [JCAMP-DX 5.01](https://doi.org/10.1351/pac199971081549)
- [EMR protocol](https://doi.org/10.1351/pac200678030613)

The source header gives the contact talos@spindynamics.org.
