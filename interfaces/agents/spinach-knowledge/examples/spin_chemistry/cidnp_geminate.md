# examples/spin_chemistry/cidnp_geminate.m

Source: [examples/spin_chemistry/cidnp_geminate.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/spin_chemistry/cidnp_geminate.m)

## Purpose

A minimal geminate CIDNP calculation that starts from an electron-pair singlet, evolves a radical-pair reaction model, and reports proton longitudinal magnetisation separately in reactants and products.

## Spin system and reaction model

The system is specified at 14.1 T with isotopes E, E, and 1H. The Zeeman scalar entries are 2.0023, 2.0024, and 1.0. The only explicit scalar coupling is between the second electron and the proton, with value 1e7. An explicit first-order loss record uses the singlet selector on electrons [1, 2], rate 1e7, and no tracked products (the triplet rate is zero). The selector gives the Haberkorn drain. The source does not annotate units for these rate and coupling entries.

The Hamiltonian is built under the ESR assumption, and the initial state is the singlet of spins 1 and 2. The script doubles the state space into reactant and product sectors, sets the product-sector Hamiltonian to zero (no product-subspace dynamics), and builds the kinetic blocks so that loss from reactants is balanced by gain in products. No relaxation superoperator is configured; the evolved Liouvillian is H + iK.

## Evolution and observable

The doubled state is evolved for one final step of 1 microsecond. The unweighted proton Lz coil is then projected onto the reactant and product portions of the state vector, and the script prints the real-valued longitudinal magnetisation for each sector. These are model outputs requested by the code; this source-only description does not assert their numerical values or interpret them as measurements.

## Citation

The source file supplies author comments but no DOI or publication citation.
