# kernel/overloads/@polyadic/allfinite.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/allfinite.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=polyadic/allfinite.m)

`allfinite(p)` visits every inner core cell in every polyadic term, then every prefix and suffix cell. It passes each stored component to `allfinite`; the first false result returns false immediately, and true is returned only after all component checks succeed. Recursive dispatch covers components that are themselves polyadic.

This is a finiteness test on the stored representation, including the factors and the left and right operators. It does not construct the Kronecker products or materialise the represented matrix. No dimension validation, matrix action, or broadcasting is performed by this method; component handling is delegated to the applicable `allfinite` overloads.