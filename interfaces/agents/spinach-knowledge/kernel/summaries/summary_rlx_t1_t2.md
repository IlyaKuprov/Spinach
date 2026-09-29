# kernel/summaries/summary_rlx_t1_t2.m

- Signature: `summary_rlx_t1_t2(spin_system,header)`

## Behaviour

Prints one extended T1/T2 row per spin, listing index, isotope, stored `r1_rates{n}` and `r2_rates{n}`, and spin label. A scalar numeric rate appears in signed scientific notation with five decimal places; a nonscalar numeric rate appears as `anisotropic`; a nonnumeric rate appears as `orientation`. The routine performs no unit conversion.

The supplied header precedes the table. Output is routed through `report` rather than returned. The local guard requires a structure `spin_system` and character-array `header`.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/summaries/summary_rlx_t1_t2.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=summary_rlx_t1_t2.m)
