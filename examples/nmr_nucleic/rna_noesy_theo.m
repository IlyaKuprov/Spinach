% 1H-1H NOESY spectrum of the example RNA molecule provided by
% the Gerhard Wagner group at Harvard University.
%
% Calculation time: hours
%
% Shunsuke Imai
% Scott Robson
% Gerhard Wagner
% Zenawi Welderufael
% Ilya Kuprov

function rna_noesy_theo()

% Import RNA data
options.noshift='delete';
options.deut_list={'GUA:H1','GUA:H21','GUA:H22','CYT:H41',...
                   'CYT:H42','URI:H3','ADE:H61','ADE:H62'};
[sys,inter]=nuclacid('example.pdb','example.txt',options);

% Magnet field
sys.magnet=17.62;

% Tolerances
sys.tols.inter_cutoff=1.0;
sys.tols.prox_cutoff=5.0;
sys.disable={'krylov','colorbar'};
sys.enable={'prop_cache','greedy'};

% Relaxation theory
inter.relaxation={'redfield'};
inter.rlx_keep='kite';
inter.equilibrium='zero';
inter.tau_c={3e-9};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.level=5; bas.space_level=3;

% Create the spin system structure
spin_system=create(sys,inter);

% Kill carbons and nitrogens (RNA assumed unlabelled)
spin_system=kill_spin(spin_system,strcmp('13C',spin_system.comp.isotopes));
spin_system=kill_spin(spin_system,strcmp('15N',spin_system.comp.isotopes));

% Build the basis
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.tmix=0.200;
parameters.offset=3473;
parameters.sweep=[7500 7500];
parameters.npoints=[512 1024];
parameters.zerofill=[1024 4096];
parameters.spins={'1H'};
parameters.axis_units='ppm';
parameters.rho0=state(spin_system,'Lz','1H','cheap');

% Simulation
fid=liquid(spin_system,@noesy,parameters,'nmr');

% Apodization
fid.cos=apodization(fid.cos,'sqcosbell-2d');
fid.sin=apodization(fid.sin,'sqcosbell-2d');

% F2 Fourier transform
f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1));
f1_sin=real(fftshift(fft(fid.sin,parameters.zerofill(2),1),1));

% States signal
f1_states=f1_cos-1i*f1_sin;

% F1 Fourier transform
spec=fftshift(fft(f1_states,parameters.zerofill(1),2),2);

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,-real(spec),parameters,...
        20,[0.01 0.5 0.01 0.5]/200,2,256,6,'positive');

end

