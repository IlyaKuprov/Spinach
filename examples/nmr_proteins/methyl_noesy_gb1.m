% 1H-1H NOESY spectrum of GB1 with everything deuterated except 
% methyl groups. Deuteria are kept in the spin system because
% they are a part of the coupling network; methyl group rotati-
% on is not accounted for in this simulation.
%
% Calculation time: hours.
%
% ilya.kuprov@weizmann.ac.il

function methyl_noesy_gb1()

% Protein data import
options.pdb_mol=1;
options.select='all';
options.noshift='delete';
options.deuterate='non-Me';
[sys,inter]=protein('2N9K.pdb','2N9K.bmrb',options);

% Magnet field
sys.magnet=21.1356;

% Tolerances
sys.tols.inter_cutoff=100;              % Only significant DD couplings
sys.tols.prox_cutoff=5.0;               % Increase till convergence

% Relaxation theory
inter.relaxation={'redfield'};
inter.rlx_keep='kite';
inter.equilibrium='zero';
inter.tau_c={5e-9};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.level=2; bas.space_level=2;         % Pairwise approximation

% Algorithmic options
sys.enable={'prop_cache','op_cache','greedy'};

% Create the spin system structure
spin_system=create(sys,inter);

% Kill carbons and nitrogens (protein assumed unlabelled)
spin_system=kill_spin(spin_system,strcmp('13C',spin_system.comp.isotopes));
spin_system=kill_spin(spin_system,strcmp('15N',spin_system.comp.isotopes));

% Build the basis
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.tmix=200e-3;
parameters.offset=750;
parameters.sweep=[3000 3000];
parameters.npoints=[512 512];
parameters.zerofill=[2048 2048];
parameters.spins={'1H'};
parameters.axis_units='ppm';
parameters.rho0=state(spin_system,'Lz','1H','cheap');

% Simulation
fid=liquid(spin_system,@noesy,parameters,'nmr');

% Apodisation
fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}});
fid.sin=apodisation(spin_system,fid.sin,{{'sqcos'},{'sqcos'}});

% F2 Fourier transform
f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1));
f1_sin=real(fftshift(fft(fid.sin,parameters.zerofill(2),1),1));

% States signal
f1_states=f1_cos-1i*f1_sin;

% F1 Fourier transform
spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2);

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,-real(spectrum),parameters,...
        20,[0.00125 0.0125 0.00125 0.0125],2,256,6,'positive');

end

