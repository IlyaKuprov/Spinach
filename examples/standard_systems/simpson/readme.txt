These input files were kindly provided by Zdenek Tosner,
Thomas Vosegaard, and Niels Chr. Nielsen (Aarhus University,
Denmark).

Use s2spinach to read one literal spinsys block: nuclei, channels,
ppm shifts, dipoles, J couplings, and quadrupolar tensors. Nonzero
shift amplitudes must end in p. Expressions and pulse programs
are not run, and the caller chooses the field and Hamiltonian
approximation.

With Spinach on the MATLAB path and this directory current:

    [sys,inter]=s2spinach('ccacb.spinsys');
    sys.magnet=9.4;                    % tesla
    bas.formalism='zeeman-hilb';
    bas.approximation={'none'};
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
