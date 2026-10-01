# kernel/overloads/@polyadic/core_specs.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/core_specs.m)

`cores=core_specs(p)` reconstructs the nested core descriptions accepted by the constructor, pairing ordinary action handles with their stored dimensions and adjoints. Numeric factors are returned unchanged. This is used internally when combining terms or passing a factor list to `kronm`; it does not execute any action.
