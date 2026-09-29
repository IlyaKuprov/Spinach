# interfaces/jsonlab-1.5/loadubjson.m

- MATLAB source: [loadubjson.m](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/jsonlab-1.5/loadubjson.m)
- Related documentation: [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=Main_Page)
- Signature: `data = loadubjson(fname,varargin)`

## Contract

Reads UBJSON bytes from a file opened in binary mode, or parses the supplied character input as UBJSON when it contains a square or curly bracket. Otherwise `fname` must identify an existing file; absent files and malformed or unsupported markers raise errors. Top-level values must be objects or arrays. This is a data conversion utility and assigns no physical units.

Objects become MATLAB structs, with names converted to valid MATLAB field names as needed. Arrays are reconstructed in MATLAB values: compatible numeric blocks can be numeric arrays, while general or heterogeneous content uses cells. The array parser handles UBJSON typed and counted-array headers, and the numeric-block decoder maps markers `i`, `U`, `I`, `l`, `L`, `d`, and `D` to `int8`, `uint8`, `int16`, `int32`, `int64`, `single`, and `double`, respectively. Boolean markers `T` and `F` become logical values; `Z` and `N` become `[]`. A single top-level cell is unwrapped before return.

## Calling forms and options

Supply an option struct or name-value pairs after `fname`, as in `loadubjson(fname,opt)` or `loadubjson(fname,'IntEndian','L')`.

- `IntEndian` defaults to `'B'`, the UBJSON big-endian integer representation; use `'L'` for little-endian integer fields.
- `SimplifyCell` defaults to `0`; `1` asks the parser to try `cell2mat`. If conversion cannot be performed, the existing representation is retained. The implementation also consults the internal `SimplifyCellArray` key when deciding whether to preserve a struct array produced by simplification; the source does not list that key among its public options.
- `NameIsString` defaults to `0`. Set to `1` for older UBJSON Specification Draft 8 or earlier data whose name tag is encoded as a string.

The source gives a round-trip example using `struct('string','value','array',[1 2 3])`, `saveubjson`, and `loadubjson`; it also demonstrates reading `examples/example1.ubj` with optional `'SimplifyCell',1`.

## Provenance

This JSONLab utility is credited to Qianqian Fang in its source, dated 2013-08-01; the source identifies the [JSONLab project](http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab) and BSD licensing. The repository's [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=Main_Page) is the linked general documentation entry. No DOI is cited in the source.
