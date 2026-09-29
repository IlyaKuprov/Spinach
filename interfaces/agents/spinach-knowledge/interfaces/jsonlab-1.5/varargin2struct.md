# interfaces/jsonlab-1.5/varargin2struct.m

[Canonical source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/jsonlab-1.5/varargin2struct.m) · [JSONLab project page](http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)

**Call:** `opt=varargin2struct(varargin)` converts the elements of a MATLAB argument cell into a scalar structure. With no inputs it returns an empty structure. It accepts alternating character-vector name/value pairs, scalar struct arguments, or a mixture of the two.

Each name/value key is lowercased before assignment; its following value is copied without type conversion. Struct arguments are merged field by field by `mergestruct`: later fields replace earlier fields with the same name, and struct field spelling is otherwise retained. For example, `varargin2struct('Tolerance',1e-6,struct('MaxIter',50))` yields fields `tolerance` and `MaxIter`.

A name must be a MATLAB character vector with a following value; a dangling name, a non-character/non-struct item, or a struct array that cannot be merged raises an error. This helper does not validate option meanings, normalise values, or infer units. The result is an options struct only; it does not apply the options to a caller.

The vendored source identifies Qianqian Fang as author and gives a BSD license (see `LICENSE_BSD.txt`).
