% 1H-1H NOESY spectrum of ubiquitin with 65 ms mixing time. It is
% assumed that the protein is not 13C- or 15N-labelled.
%
% Calculation time: hours.
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk

function noesy_ubiquitin()

% Protein data import
options.pdb_mol=1;
options.select='all';
options.noshift='delete';
[sys,inter]=protein('1D3Z.pdb','1D3Z.bmrb',options);

% Magnet field
sys.magnet=21.1356;

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

% Algorithmic options
sys.enable={'prop_cache','greedy'};

% Create the spin system structure
spin_system=create(sys,inter);

% Kill carbons and nitrogens (protein assumed unlabelled)
spin_system=kill_spin(spin_system,strcmp('13C',spin_system.comp.isotopes));
spin_system=kill_spin(spin_system,strcmp('15N',spin_system.comp.isotopes));

% Build the basis
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.tmix=0.065;
parameters.offset=4250;
parameters.sweep=[11750 11750];
parameters.npoints=[512 512];
parameters.zerofill=[2048 2048];
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
spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2);

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,-real(spectrum),parameters,...
        20,[0.00125 0.0125 0.00125 0.0125],2,256,6,'positive');

end

