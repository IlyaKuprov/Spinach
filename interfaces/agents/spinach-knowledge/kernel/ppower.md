# kernel/ppower.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/ppower.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ppower.m)

- Signature: `P=ppower(spin_system,P,N)`

## Behaviour

Raises a finite numeric square matrix `P` to the non-negative integer power `N`. The source requires `spin_system.tols.prop_chop`; it rejects powers above `flintmax` and converts the accepted power to `uint64` for its bit loop.

- `N=0`: returns a dense or sparse identity matching the input matrix's storage class.
- `N=1`: returns the validated input matrix unchanged.
- Larger powers: binary exponentiation multiplies selected powers into an identity accumulator and squares the working matrix as needed. Each multiplication is passed through `clean_up` using `spin_system.tols.prop_chop`.

The function has no pulse phase, amplitude, timing, or file-output interface; it returns the powered matrix.
