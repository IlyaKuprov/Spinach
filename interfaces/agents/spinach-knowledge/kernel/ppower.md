# kernel/ppower.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/ppower.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ppower.m)

`P=ppower(spin_system,P,N)` raises a finite numeric square matrix to a non-negative integer power. Powers above `flintmax` are rejected. For `N=0`, the result is an identity matching the input dimensions and sparse/full representation.

Binary exponentiation uses panel multiplication with `spin_system.tols.prop_chop`. Ordinary cleanup owns rounding and storage decisions; results can differ from a single native matrix power by numerical roundoff. Panel processing targets transient product fill-in, not a hard bound on total memory.
