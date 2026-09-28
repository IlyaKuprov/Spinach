# interfaces/jsonlab-1.5/jsonopt.m

- Signature: `val=jsonopt(key,default,varargin)`

Looks up an option in the first optional argument, expected to be a struct. It returns the field matching `key` exactly; if absent, it tries `lower(key)`. If neither field exists, the argument is not a struct, or no optional argument is supplied, it returns `default`. Additional optional arguments are not used.

The option struct can be produced by `varargin2struct` from parameter/value pairs. Part of the JSONLab toolbox ([project page](http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)).
