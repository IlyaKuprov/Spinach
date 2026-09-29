# interfaces/jsonlab-1.5/jsonopt.m

[Canonical source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/jsonlab-1.5/jsonopt.m) · [JSONLab project page](http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)

**Call:** `val=jsonopt(key,default,varargin)`.

The function starts with `val=default`; if no optional argument is supplied, it returns that default. Otherwise it examines only the first value in `varargin`. If that value is a struct, it first checks for a field named exactly `key`; if absent, it checks for a field named `lower(key)`. It returns the matching field value unchanged, or `default` when the optional value is not a struct or neither field exists. Further optional arguments are ignored.

The JSONLab source notes that an options struct can be built with `varargin2struct` from parameter/value pairs. This helper does not merge options, coerce field values, or validate the requested key; field-name lookup follows MATLAB `isfield` semantics.
