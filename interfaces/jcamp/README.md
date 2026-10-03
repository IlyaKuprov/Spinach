# JCAMP export of NMR and EPR data

```matlab
text=jcamp_export(data);
```

One structure describes the entire file. The result is an ASCII character row;
if `data.filename` is present, the same text is also written to that file.
There is no JSONLab dependency. Existing files are replaced only after the
complete output has been serialised and written successfully. A directory
cannot be used as `data.filename`.

## Export directly from Spinach results

The wrappers build `data` internally and call `jcamp_export` as the final
writer. Use them with the same `spin_system`, `parameters`, and result arrays
as the example calculations; no manual variables or pages are needed.

```matlab
info.title='Proton pulse acquisition';
info.origin='Your institution';
info.owner='Your name';
info.filename='proton.jdx';    % '' returns text only
info.sequence='Pulse acquisition';
info.metadata=cell(0,2);      % additional JCAMP records, if needed
info.delay=[0 0];            % simulated RD and ID in microseconds
info.acquisition='SIMULTANEOUS';
text=jcamp_nmr(spin_system,parameters,fid,{'time'},info);
text=jcamp_nmr(spin_system,parameters,spectrum,{'frequency'},info);
```

The two delay values and acquisition convention are explicit scientific
metadata, not inferred from a sequence name. Frequency-only results do not
require them. The time axes start at acquisition zero; pre-acquisition evolution
such as `parameters.dead_time` belongs in the delay/sequence description.
`info.metadata` can add sample, referencing, and processing records; generated
records cannot be overridden. Use ASCII descriptions, as required by JCAMP.

### NMR arrays and quadrature structures

`jcamp_nmr(spin_system,parameters,signal,domains,info)` supports:

- 1D column FIDs and spectra, as in `acquire` and `plot_1d`.
- 2D `[F2,F1]` arrays, as in COSY and `plot_2d`; `fid.pos/fid.neg` from
  HSQC and `fid.cos/fid.sin` from NOESY are accepted directly.
- 3D `[F1,F2,F3]` arrays, as in HNCO, HNCACO, and `plot_3d`; the four
  `pos_pos/pos_neg/neg_pos/neg_neg` components are accepted directly.
- Mixed-domain arrays and singleton time dimensions. `domains` always lists
  the **physical** F1/F2/F3 order, even for the reversed 2D array order.

```matlab
text=jcamp_nmr(spin_system,parameters,fid,{'time','time'},info);
text=jcamp_nmr(spin_system,parameters,spectrum,{'frequency','frequency'},info);
text=jcamp_nmr(spin_system,parameters,fid,{'time','time','time'},info);
text=jcamp_nmr(spin_system,parameters,spectrum,...
               {'frequency','frequency','frequency'},info);
```

`parameters.sweep`, `offset`, and `spins` follow physical-dimension order.
A scalar entry is shared across dimensions, matching homonuclear examples.
Array dimensions determine sample counts; no resizing to `npoints` or
`zerofill` occurs. All component arrays must have the same declared shape.
Named components become separate LINK blocks with their names in titles and
`$SPINACH COMPONENT`; complex values become signed real/imaginary pages.
No echo/antiecho or States recombination occurs: those operations are
sequence-specific and remain in the caller's processing code.

Time coordinates are `(0:N-1)/sweep`. Frequency coordinates use `ft_axis`,
including odd/even FFT lengths and the non-duplicated periodic edge, exactly
as Spinach's spectral plots. These frequency dimensions require at least
three points, following `ft_axis`'s contract. Exported frequency axes are Hz,
not `parameters.axis_units` display conversions; observation frequencies are
computed from the field and nuclei in MHz. Private axis records identify
all nuclei, frequencies, offsets, and domains in physical order. The ordinary
observe record identifies the directly detected nucleus. Mixed-domain blocks
use the data type of their tabulated, first-array-dimension axis.

### EPR acquisitions, spectra, and scans

`jcamp_epr(spin_system,parameters,signal,kind,info)` covers common outputs.
EMR uses the same ownership/file fields, plus explicit method metadata:

```matlab
info.title='Nitroxide field sweep';
info.origin='Your institution';
info.owner='Your name';
info.filename='nitroxide.jdx';
info.metadata=cell(0,2);
info.detection='CW';
info.method='SPECTRUM';
info.description='Nitroxide model and numerical settings used for this run';
[spec,parameters]=fieldsweep(spin_system,parameters);
text=jcamp_epr(spin_system,parameters,spec,'field',info);
```

Kinds and native shapes:

- `'field'`: `fieldsweep` row result and its **returned** `parameters.b_axis`
  row in tesla; `parameters.mw_freq` is recorded in Hz.
- `'endor'`: `endor_davies`/`endor_mims` row result and `parameters.n_frq`
  row in Hz. Signed RF coordinates are preserved, not absolute-value folded.
  Set `info.detection='PULSE'` and `info.method='ENDOR'`.
- `'time'`: column pulse-acquire FID or `[F2,F1]` pulse array; dimensionality
  comes from `parameters.npoints`, or `parameters.nsteps` for HYSCORE.
  Dwell times are `1/parameters.sweep`. Set the actual method explicitly.
- `'frequency'`: already processed column spectrum or `[F2,F1]` spectrum;
  dimensionality comes from `parameters.zerofill`, actual lengths from the
  signal, and axes from `sweep`, `offset`, and `ft_axis`.

```matlab
text=jcamp_epr(spin_system,parameters,answer,'endor',info);
text=jcamp_epr(spin_system,parameters,fid,'time',info);
text=jcamp_epr(spin_system,parameters,spectrum,'frequency',info);
```

The scan wrappers explicitly translate the native row storage into column
traces; they reject other shapes. Regular acquisition wrappers do not transpose
or reshape their input. Named component structures also work. Results are
`EMR SIMULATION`, with Spinach as the simulation source and the supplied
physical description as simulation parameters. Method names are JCAMP core
identifiers, not MATLAB sequence names. For example, DEER uses `ELDOR` and
its particular pulse sequence is described in the metadata. The wrappers
never infer detection mode, pulse timings, reference conventions, or
instrument settings merely from an array.

### Explicit sampling and other magnetic resonance data

For DEER/ESEEM delay traces, echo trajectories, irregular time grids, or
custom EMR maps, use physical axes directly:

```matlab
info.detection='PULSE';
info.method='ELDOR';
info.description='Four-pulse DEER; actual model and pulse settings';
text=jcamp_signal(spin_system,{delays_seconds},deer_signal,...
                  {'SECOND'},{'DEER delay'},info);
```

`jcamp_signal(spin_system,axes,signal,units,names,info)` accepts one to three
axis **columns**, in MATLAB array-dimension order. Units and names are cell
rows with one entry per axis. A 1D signal is a column; in higher dimensions,
`size(signal,n)==numel(axes{n})`. Explicit EMR axes use `SECOND`, `HERTZ`,
`TESLA`, and the other protocol keywords, not plot display units such as
microseconds or MHz. Numeric arrays and named component structures are both
accepted. This covers custom ENDOR/HYSCORE/ELDOR/ESEEM, angular scans, and
other sampled EMR results within the EMR unit/method vocabulary, without
inventing a new sampling convention. Diagnostic routines that only plot and
return no results cannot be exported until the caller obtains their arrays.

For non-standard NMR sampling/layouts, `jcamp_grid(axes,signal,units,names,block)`
is the shared sampled-array translator: supply a typed NMR block with its
metadata, then put the returned blocks into `jcamp_export`'s file structure.
Use the low-level writer directly for peak tables, assignments, or ragged
pages. A JCAMP magnetic-resonance file is not a general replacement for MRI
image formats or a format for spin-state trajectories with no observable.

All wrapper outputs share the precision and reader-compatibility limitations
below. Multidimensional JCAMP support is reader-dependent; preserving all
Spinach components does not imply a vendor can reconstruct their quadratures.
MATLAB `jcampread` recovered the native 1D FID/spectrum in the wrapper checks,
but rejected the linked NOESY output. Not all native explicit-pair or
multidimensional NTUPLES were recovered by nmrglue and jcampconverter either.
Use a reader whose support covers the required layout; these wrappers do not
claim universal vendor import compatibility.

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
  A period may prefix a technique-specific label; embedded periods are
  permitted only in the user-defined private namespace.
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
retain their explicit units. NMR ordinate units are `ARBITRARY UNITS`,
`MAGNITUDE`, or `POWER` (JCAMP 5.01); EMR ordinate units are `ARBITRARY UNITS`,
`INTENSITY`, or `POWER`. These units are also checked for each dependent
page variable. Complex NMR traces require `ARBITRARY UNITS`; magnitude/power
labels describe already-real transformed data. The exporter never computes a magnitude or power merely
because the corresponding label is supplied, and performs no unit conversion
or processing.

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
Pages retain their input order, including repeated coordinate sets. To label
replicates for readers keyed by coordinate, declare an additional `PAGE`
variable and fix its value in every page.
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
- `width`: optional non-negative finite column, in `xunits`; unassigned and
  EMR widths require `y`;
- `multiplicity`: optional NMR-only cell column of `S`, `D`, `T`, `Q`, `M`, or `U`,
  requiring `y` only in unassigned tables;
- `assignment`: optional cell column of ASCII assignment strings;
- `method`: ASCII description of peak finding and width convention, required
  whenever width or EMR assignment is supplied. NMR assignments without
  widths do not require a method comment.

Without assignments, a `PEAK TABLE` contains heights and optionally either
width or multiplicity. With assignments, `PEAK ASSIGNMENTS` contains angle-
bracketed strings inside parenthesised groups. NMR uses `DATA CLASS=ASSIGNMENTS`
and the `PEAK ASSIGNMENTS` table label; EMR uses `DATA CLASS=PEAK ASSIGNMENTS`
as defined in Section 4.1.4 (the protocol's summary table instead lists
`ASSIGNMENTS`).
Unassigned EMR lists are `(XY)` or `(XYW)`; NMR lists use repeated markers,
such as `(XY..XY)` or `(XYW..XYW)`.
NMR widths and multiplicities may both be supplied, in the standard's `XYMWA`
order. Both NMR and EMR permit x-only assignments. NMR assigned heights, widths,
and multiplicities are independently optional; EMR widths require heights. Assignment strings cannot contain
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
`.SIMULATION PARAMETERS` for simulations; both descriptions must be non-empty
ASCII text, either a row string or a cell vector of lines, not numeric arrays.
Acquisition mode is `SIMULTANEOUS`,
`SEQUENTIAL`, or `SINGLE`; EMR detection is `CW` or `PULSE`. The EMR method
is one of the protocol's core identifiers: `DYNAMIC`, `ELDOR`, `ENDOR`, `ESEEM`,
`ODMR`, `GONIOMETER`, `HYSCORE`, `KINETIC`, `SATURATION`, `SPECTRUM`, `FID`,
`TRIPLE`, `IMAGING`, or `SPECTRAL SPATIAL`. Additional experiment abbreviations
can be described using referenced private metadata as suggested by the EMR
protocol, while selecting the applicable core method.

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
