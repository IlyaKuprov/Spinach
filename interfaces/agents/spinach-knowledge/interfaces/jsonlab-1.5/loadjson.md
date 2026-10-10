# interfaces/jsonlab-1.5/loadjson.m

- MATLAB source: [loadjson.m](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/jsonlab-1.5/loadjson.m)
- Related documentation: [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=Main_Page)
- Signature: `data = loadjson(fname,varargin)`

## Contract

Reads JSON text from a file or from the supplied character input and converts it to MATLAB values. The input is treated as inline text when the source's JSON-looking pattern matches; otherwise `fname` must name an existing file. File contents are read with `fileread` (with a `file://` fallback). Missing files and malformed syntax raise errors. The parser expects each top-level value to begin with an object or array delimiter.

JSON objects become MATLAB structs; property names are converted to valid MATLAB field names when needed. Arrays are parsed as MATLAB arrays where the contents can be combined, or as cells when their elements do not combine. JSON strings become character data, numbers become numeric values, booleans become logical scalars, and `null` becomes `[]`. A single top-level parsed cell is unwrapped before return. This is a data conversion utility, not a physical measurement interface: it assigns no units.

## Calling forms and options

Pass an option struct or name-value pairs after `fname`, for example `loadjson(fname,opt)` or `loadjson(fname,'SimplifyCell',1)`.

- `SimplifyCell` defaults to `0`; when set to `1`, the parser tries `cell2mat` on parsed cell arrays. Conversion is guarded by `try/catch`, so incompatible data remain unsimplified.
- `FastArrayParser` defaults to `1`; `0` selects the legacy array parser. Values above one set the nested-array depth threshold for the fast path. The source documents cell/array shape outcomes for nested arrays; do not assume every irregular JSON array becomes a rectangular numeric matrix.
- `ShowProgress` defaults to `0`; setting it to `1` opens a MATLAB progress bar while parsing.

The implementation also reads internal `JSONLAB_ArrayDepth_` and `SimplifyCellArray` keys during recursive parsing and cell simplification. They are not listed in the source public option block and are not described here as stable caller options.

The source example demonstrates nested object access and a numeric array:

`loadjson('{"obj":{"string":"value","array":[1,2,3]}}')`

It also shows loading `examples/example1.json`, optionally with `'SimplifyCell',1`.

## Provenance

This is the JSONLab utility by Qianqian Fang, with earlier contributions credited in its source to [Nedialko Krouchev](http://www.mathworks.com/matlabcentral/fileexchange/25713), [Francois Glineur](http://www.mathworks.com/matlabcentral/fileexchange/23393), and Joel Feenstra. The source identifies the JSONLab project at [iso2mesh / JSONLab](http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab) and states a BSD license. No DOI is cited in the source.
