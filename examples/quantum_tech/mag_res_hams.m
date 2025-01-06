% Drops a few common Magnetic Resonance evolution generators, initial 
% states, and detection states. Generators are split into Zeeman, cou-
% pling, and dissipation parts.
%
% ilya.kuprov@weizmann.ac.il

function mag_res_hams() %#ok<*NASGU>

% Get strychnine data
[sys,inter]=strychnine({'1H'});

% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=4;

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={200e-12};

% Proximity cut-off
sys.tols.prox_cutoff=4.0;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Interaction representation Hamiltonian
H_full=hamiltonian(assume(spin_system,'nmr'));
H_field_part=hamiltonian(assume(spin_system,'nmr','zeeman'));
H_inter_part=hamiltonian(assume(spin_system,'nmr','couplings'));

% Relaxation superoperator
R_full=relaxation(spin_system); 

% Initial and detection state
rho=state(spin_system,'Lx','1H');
coil=state(spin_system,'Lx','1H');

% Save to disk and clear workspace
save('strychnine.mat','-v7.3'); clear('variables');

% ============================================================

% Get strychnine data
[sys,inter]=allyl_pyruvate({'1H'});

% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=4;

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={200e-12};

% Proximity cut-off
sys.tols.prox_cutoff=4.0;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Interaction representation Hamiltonian
H_full=hamiltonian(assume(spin_system,'nmr'));
H_field_part=hamiltonian(assume(spin_system,'nmr','zeeman'));
H_inter_part=hamiltonian(assume(spin_system,'nmr','couplings'));

% Relaxation superoperator
R_full=relaxation(spin_system); 

% Initial and detection state
rho=state(spin_system,'Lx','1H');
coil=state(spin_system,'Lx','1H');

% Save to disk and clear workspace
save('allyl_pyruvate.mat','-v7.3'); clear('variables');

% ============================================================

% Protein data import
options.pdb_mol=1;
options.select='all';
options.noshift='delete';
[sys,inter]=protein('..\nmr_proteins\1D3Z.pdb',...
                    '..\nmr_proteins\1D3Z.bmrb',options);

% Magnet field
sys.magnet=14.1;

% Tolerances
sys.tols.inter_cutoff=2.0;
sys.tols.prox_cutoff=4.0;

% Relaxation theory
inter.relaxation={'redfield'};
inter.rlx_keep='kite';
inter.equilibrium='zero';
inter.tau_c={5e-9};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.level=4; bas.space_level=3;

% Create the spin system structure
spin_system=create(sys,inter);

% Kill carbons and nitrogens (protein assumed unlabelled)
spin_system=kill_spin(spin_system,strcmp('13C',spin_system.comp.isotopes));
spin_system=kill_spin(spin_system,strcmp('15N',spin_system.comp.isotopes));

% Build the basis
spin_system=basis(spin_system,bas);

% Interaction representation Hamiltonian
H_full=hamiltonian(assume(spin_system,'nmr'));
H_field_part=hamiltonian(assume(spin_system,'nmr','zeeman'));
H_inter_part=hamiltonian(assume(spin_system,'nmr','couplings'));

% Relaxation superoperator
R_full=relaxation(spin_system); 

% Initial and detection state
rho=state(spin_system,'Lx','1H');
coil=state(spin_system,'Lx','1H');

% Save to disk 
save('ubiquitin.mat','-v7.3');

end

