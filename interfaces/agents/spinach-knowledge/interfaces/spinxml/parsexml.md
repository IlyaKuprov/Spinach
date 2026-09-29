# interfaces/spinxml/parsexml.m

[MATLAB implementation](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/spinxml/parsexml.m) · [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=parsexml.m)

## Interface and result structure

`xml=parsexml(filename)` reads a named XML file with MATLAB's `xmlread(filename,'XMLEngine','maxp')` and recursively converts the resulting DOM children to MATLAB structures. The implementation is used by `x2spinach()` for SpinXML imports; its source comment discourages direct calls. `filename` must be a non-empty character array naming an existing file. Read failures and conversion failures are rethrown as errors naming the file.

The result is a structure array with one entry per child node, in child order; an empty child list is represented by `[]`. Each entry has four fields: `name`, `attributes`, `data`, and `children`. Element names come from `TagName`; text and comment nodes are named `#text` and `#comment`, and other nodes use their MATLAB DOM class name. Element-node `data` is empty; nodes with `TextContent` store it as a character array, and otherwise `data` is empty. `children` recursively has the same representation.

`attributes` is `[]` when the node has none. Otherwise it is a structure array whose entries each have `name` and `value` character fields, copied from the DOM attributes. These are structural XML values: the parser defines no Spinach state/operator dimensions, unit conversion, or caching behaviour.
