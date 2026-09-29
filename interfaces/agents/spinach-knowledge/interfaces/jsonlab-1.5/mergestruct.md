# interfaces/jsonlab-1.5/mergestruct.m

- MATLAB source: [mergestruct.m](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/jsonlab-1.5/mergestruct.m)
- Related documentation: [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=Main_Page)
- Signature: `s = mergestruct(s1,s2)`

## Contract

Combines the fields of two MATLAB struct inputs into one struct. Both inputs must satisfy `isstruct`; inputs with `length>1` are rejected, so struct arrays are not accepted. The implementation starts from `s1` and assigns every field from `s2` in turn: a same-named field in `s2` replaces the value in `s1`, and fields present only in `s2` are added. This is a shallow field-wise replacement, not a recursive merge; nested structs are copied as field values. Empty struct inputs pass the explicit array-length guard. No units or numerical transformations are applied.

## Example behaviour

For `s1 = struct('a',1,'b',2)` and `s2 = struct('b',3,'c',4)`, the returned struct has `a = 1`, `b = 3`, and `c = 4`. The implementation raises an error for non-struct inputs and for struct inputs whose length exceeds one.

## Provenance

The source credits Qianqian Fang and dates the function 2012-12-22. It identifies the [JSONLab project](http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab) and BSD licensing. The repository's [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=Main_Page) is the linked general documentation entry. No DOI is cited in the source.
