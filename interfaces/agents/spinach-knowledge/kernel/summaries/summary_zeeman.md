# kernel/summaries/summary_zeeman.m

- Signature: `summary_zeeman(spin_system,header)`

## Behaviour

For each spin with a nonempty `spin_system.inter.zeeman.matrix{n}`, prints its index, isotope, multiplicity, three matrix rows, isotropic value `trace(M)/3`, and spectral matrix norms of the rank-1 and rank-2 parts. The parts come from `mat2sphten(M)` and `sphten2mat`; each norm uses MATLAB `norm(...,2)`. Values are printed as stored, without a unit conversion in this routine.

The supplied `header` precedes the report. There is no return value; output goes through `report`. The local guard requires a structure `spin_system` and character-array `header`.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/summaries/summary_zeeman.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=summary_zeeman.m)
