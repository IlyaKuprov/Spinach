# fme_symbolic_check.m

- Source: `examples/fundamentals/fme_symbolic_check.m`
- Entry point: `fme_symbolic_check()`; no inputs or outputs
- Requires Symbolic Math Toolbox
- Reference: Bock, Ernst, and Uhrig, DOI:10.5194/mr-2026-10, Appendix A and supplement Section S5

Spinach pauli(2) supplies the spin operators, which are converted to exact symbolic matrices. The resulting noncommuting spin-half model checks the mean Hamiltonian and first inverse-frequency AHT correction by evaluating the first two Magnus integrals. Independent Fourier commutators are required to agree symbolically. Integration of the oscillatory Hamiltonian supplies the leading AHT kick; removal of its mean supplies the Van Vleck kick. Their Fourier coefficients are obtained by exact integrals, rather than fitting values at selected phases.

After substituting a=3, b=2, c=5 and omega=17, the symbolic effective generator and every leading kick coefficient are compared with floquet_magnus. Symbolic identities must vanish exactly; numerical differences must lie within a documented floating-point bound. The example validates low orders for this model, not all higher orders or the supplement's MAS calculations. No plot is produced.
