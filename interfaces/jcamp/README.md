# JCAMP export of NMR and EPR data

```matlab
text=jcamp_export(data);
```

One structure describes the entire file. The result is an ASCII character row;
if `data.filename` is present, the same text is also written to that file.
There is no JSONLab dependency. Existing files are replaced only after the
complete output has been serialised and written successfully. A directory
cannot be used as `data.filename`.

## File and block structure

Required file fields are `title`, `origin`, `owner` (non-empty ASCII character
rows), and `blocks` (a non-empty cell vector of scalar structures).
`title` identifies a compound file; each data block has its own required
`title`. A single block is written directly; two or more are enclosed in a
JCAMP `LINK` block, with a block count and unique positive `BLOCK_ID` values.
`origin` and `owner` apply to every block and are never guessed.

Every block also has:

- `type`: `NMR FID`, `NMR SPECTRUM`, `NMR PEAK TABLE`,
  `NMR PEAK ASSIGNMENTS`, `EMR SIMULATION`, or `EMR MEASUREMENT`.
  The EMR standard covers EPR/ESR, ENDOR, ESEEM, HYSCORE, and other methods;
  it does **not** use `EPR SPECTRUM` as the data type.
- `metadata`: an N-by-2 cell array of label/value pairs, without `##` or `=`.
  Use `cell(0,2)` for an empty table, though technique-specific identifiers
  described below are required. A value is a printable ASCII character row,
  a cell vector of ASCII lines, or a finite real numeric vector. Numbers
  are written as AFFN; numeric NMR `.DELAY` pairs acquire parentheses.
  Textual NMR delays must be a finite real numeric pair `(RD, ID)`.
  Generic, technique-specific (`.`), and private (`$`) labels are preserved.
  Duplicate labels, including aliases differing only in spaces, dashes,
  slashes, underscores, or case, are refused. Generated structural labels
  and numeric scaling attributes cannot be overridden.

Supply exactly one of the following data representations in each block.

## A trace: `x`, `y`, `xunits`, `yunits`

`x` and `y` are non-empty, equally sized floating-point **column vectors**.
`x` is finite and real; `y` may be real or complex and may contain NaN for
missing observations. Infinity is refused. Units are explicit non-empty ASCII strings without commas.
Optional `xname` and `yname` label the axes.

An axis exactly matching `linspace(x(1),x(end),numel(x))'` with a finite first
ordinate uses `XYDATA`; otherwise, including a singleton or leading missing
ordinate, it uses `XYPOINTS`. This conservative choice
avoids fitting or rounding genuinely irregular coordinates. Ascending,
descending, and non-monotonic explicit coordinates retain their input order.
Complex traces use NTUPLES with separate `R` and `I` pages and a page counter
`N`. Real and imaginary signs are preserved; no magnitude calculation occurs.

NMR frequencies use `HZ`, FID times use `SECONDS`, and peak positions may use
`PPM` or `HZ` under the NMR protocol. Store a spectrum in Hz and supply
`.SHIFT REFERENCE` and `.OBSERVE FREQUENCY` for chemical-shift display rather
than silently changing the abscissa units. EMR axis keywords include `TESLA`,
`HERTZ`, `SECOND`, `DEGREE`, `KELVIN`, and `WATT`. EMR microwave frequencies
are in **Hz**, unlike NMR observation frequencies, which are in **MHz**.
Tabulated abscissa units are checked against the declared data type, including
the abscissa variable used by each NTUPLES page. Fixed coordinate variables
retain their explicit units. The exporter performs no unit conversion or processing.

Example: a complex FID already calculated by Spinach:

```matlab
block=struct();
block.title='Proton pulse-acquire';
block.type='NMR FID';
block.metadata={'.OBSERVE FREQUENCY',observe_mhz;
                '.OBSERVE NUCLEUS','^1H';
                '.DELAY',[0 0];
                '.ACQUISITION MODE','SIMULTANEOUS';
                '.PULSE SEQUENCE','Pulse Acquisition'};
block.x=(0:numel(fid)-1)'/sweep_hz;
block.y=fid;
block.xunits='SECONDS';
block.yunits='ARBITRARY UNITS';
data.title=block.title;
data.origin='Your institution';
data.owner='Your name';
data.blocks={block};
data.filename='proton_fid.jdx';
text=jcamp_export(data);
```

The zero delays describe this simulated acquisition; measured data need the
actual pre-acquisition delays in microseconds. `fid` must already be a column.
No signal reshaping, Fourier transform, apodisation, normalisation, quadrature
reconstruction, or referencing is performed by the exporter.

## General NTUPLES: `variables` and `pages`

This representation is intended for multidimensional, hypercomplex,
multichannel, mixed-domain, and irregularly sampled data. It also accommodates
ragged pages with different point counts and abscissae. Dimension order and
quadrature conventions are explicit, rather than inferred from a MATLAB array.

`variables` is a non-empty structure vector. Every element has:

- `name`: non-empty ASCII text without commas;
- `symbol`: a unique uppercase identifier starting with a letter, followed
  by letters or digits, e.g. `F1`, `F2`, `R`, `I`, `N`;
- `type`: `INDEPENDENT`, `DEPENDENT`, or `PAGE`;
- `units`: ASCII text without commas; a dimensionless page counter may use `''`.

`pages` is a non-empty structure vector. Every element has:

- `x`, `y`: equally sized real floating-point column vectors; finite `x`,
  finite or NaN `y`;
- `xvar`, `yvar`: symbols naming an independent and a dependent variable;
- `coordinates`: an N-by-2 cell array of symbol/finite-scalar pairs fixing
  **every** remaining independent or page variable, with no duplicates.
  Integer coordinates must lie within the exact double-integer range
  (`abs(value)<=flintmax`); floating-point coordinates use their stored values.

A page table contains one ordinate component. Supply separate pages and
separately named dependent variables for real/imaginary, cosine/sine,
echo/antiecho, receiver channels, etc. All declared variables must occur.
Attribute lists (`VAR_DIM`, `FIRST`, `LAST`, `MIN`, `MAX`, `FACTOR`) are derived
from the supplied samples and coordinates. A tabulated variable's `VAR_DIM`
is its maximum page length; a coordinate variable's is the number of distinct
coordinates. Every page also carries its actual `NPOINTS`. The table is
incremental only for an exactly regular axis whose length and endpoints agree
with the shared variable attributes and a finite initial ordinate; otherwise
explicit pairs with the NTUPLES `PROFILE` display method are used. This
prevents shared NTUPLES attributes from changing a page-specific sampling grid.

Example: a real 2D spectrum `spectrum` with rows indexed by `f2_hz` and columns
indexed by `f1_hz`:

```matlab
block=struct();
block.title='Two-dimensional spectrum';
block.type='NMR SPECTRUM';
block.metadata={'.OBSERVE FREQUENCY',observe_mhz;
                '.OBSERVE NUCLEUS','^1H';
                '.PULSE SEQUENCE',sequence_description};
block.variables=struct('name',{'Indirect frequency','Direct frequency','Intensity'},...
                       'symbol',{'F1','F2','Y'},...
                       'type',{'INDEPENDENT','INDEPENDENT','DEPENDENT'},...
                       'units',{'HZ','HZ','ARBITRARY UNITS'});
for n=1:numel(f1_hz)
    block.pages(n).x=f2_hz;
    block.pages(n).y=spectrum(:,n);
    block.pages(n).xvar='F2';
    block.pages(n).yvar='Y';
    block.pages(n).coordinates={'F1',f1_hz(n)};
end
```

Add further independent coordinate variables for further dimensions. For
complex pages, declare `R` and `I`, and give each component a page at each
indirect coordinate. For a multichannel 1D dataset, declare a `PAGE` variable
`N` and supply its value in every page. Describe the acquisition/processing
conventions in the metadata. Generic NTUPLES can encode these data, but the
1999 specification explicitly notes the absence of a validated multidimensional
NMR profile; individual readers may support only a subset of these layouts.

## Peaks: `peaks`, `xunits`, `yunits`

`peaks` is a scalar structure with `x` (finite real floating-point column) and:

- `y`: finite real peak heights, same column shape;
- `width`: optional non-negative finite column, in `xunits`, requiring `y`;
- `multiplicity`: optional NMR-only cell column of `S`, `D`, `T`, `Q`, `M`, or `U`,
  requiring `y`;
- `assignment`: optional cell column of ASCII assignment strings;
- `method`: ASCII description of peak finding and width convention, required
  whenever width or assignment is supplied.

Without assignments, a `PEAK TABLE` contains heights and optionally either
width or multiplicity. With assignments, `PEAK ASSIGNMENTS` contains angle-
bracketed strings inside parenthesised groups. NMR uses `DATA CLASS=ASSIGNMENTS`
and the `PEAK ASSIGNMENTS` table label; EMR uses `DATA CLASS=PEAK ASSIGNMENTS`
as defined in Section 4.1.4 (the protocol's summary table instead lists
`ASSIGNMENTS`).
Unassigned EMR lists are `(XY)` or `(XYW)`; NMR lists use repeated markers,
such as `(XY..XY)` or `(XYW..XYW)`.
NMR widths and multiplicities may both be supplied, in the standard's `XYMWA`
order. NMR assignments require
heights. EMR also permits x-only assignments. Assignment strings cannot contain
angle brackets; parentheses and commas are permitted inside the brackets. Atom-number assignments require the appropriate
`CROSS REFERENCE` metadata pointing to a separately available chemical structure;
this exporter does not create molecular structures or invent assignments.

The presence of `assignment` must agree with the NMR peak data type. A peak
finding/width convention is emitted as a `$$` comment immediately before the
peak rows, so it cannot terminate the peak data record.

## Technique metadata

The exporter checks the minimal technique identifiers: NMR observation
frequency and nucleus; additionally `.DELAY` and `.ACQUISITION MODE` for FIDs;
EMR `.DETECTION MODE` and `.METHOD`; additionally `.SIMULATION SOURCE` and
`.SIMULATION PARAMETERS` for simulations. Acquisition mode is `SIMULTANEOUS`,
`SEQUENTIAL`, or `SINGLE`; EMR detection is `CW` or `PULSE`.

This is not an experiment-completeness validator. Supply the full metadata
required by the applicable IUPAC method, including where applicable:

- NMR solvent/reference records, transmitter offsets, pulse sequences,
  acquisition delays, and hypercomplex acquisition conventions.
- EMR detection method, microwave frequency/power/phase, receiver gain, scan
  time, and number of scans; CW modulation units/amplitude/frequency,
  receiver harmonic, and detection phase; ELDOR's second microwave source;
  ENDOR's static field and scanned RF power; TRIPLE's RF sources; imaging
  gradients; goniometer angle; simulation source and physical parameters.

Record sample description, instrument, long date, audit trail, processing,
references, and private labels in `metadata`. Multiline text is preserved as
continuation records, with wrapping at spaces. Reserved marker sequences `##`
and `$$` are refused inside caller-supplied text to prevent record injection.
A token that cannot fit in an 80-character line is refused, not truncated.

## Numerical encoding and compatibility

Integer metadata is written in full decimal precision; coordinates whose
integer values cannot be represented exactly as doubles are refused.
Output uses JCAMP-DX 5.01, printable ASCII, CRLF line endings, and at most
80 characters per line. AFFN uses 17 significant digits, sufficient for IEEE
binary64 round trips, and unit scale factors. Integer-only ASDF compression is
not used: no quantisation, scaling loss, or third-party codec is introduced.
Missing observations use `?` in the tables; unavailable scalar statistics are
omitted and NTUPLES attribute entries are empty. Some readers do not handle
missing observations, irregular NTUPLES, or multidimensional data.

The output is buffered in memory before writing. It is intentionally an
exporter, not a JCAMP importer or a spectrometer-vendor compatibility shim.

## Standards and reference implementations

- [Davies and Lampen, JCAMP-DX for NMR (1993)](https://iupac.org/wp-content/uploads/2021/08/JCAMP-DX_NMR_1993.pdf)
- [Lampen et al., JCAMP-DX 5.01 (1999)](https://doi.org/10.1351/pac199971081549)
- [Cammack et al., JCAMP-DX for EMR (2006)](https://doi.org/10.1351/pac200678030613)
- [McDonald and Wilks, base JCAMP-DX protocol (1988)](https://iupac.org/wp-content/uploads/2021/08/JCAMP-DX_IR_1988.pdf)
- [Lampen et al., generic NTUPLES display methods in the MS protocol (1994)](https://iupac.org/wp-content/uploads/2021/08/JCAMP-DX_MS_1994.pdf)
- [IUPAC's original reference files](https://github.com/IUPAC/JCAMP-DX)
- [nmrglue's NMR JCAMP reader](https://github.com/jjhelmus/nmrglue/blob/master/nmrglue/fileio/jcampdx.py)
- [nzhagen's JCAMP reader/writer](https://github.com/nzhagen/jcamp)
- [jcampconverter's original implementation](https://github.com/lpatiny/cheminfo-jcampconverter)

The function is written independently; reference implementation code is not
vendored into Spinach.
