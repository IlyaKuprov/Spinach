% 3D HCCH TOCSY experiment on GB1 protein.
%
% Calculation time: hours.
%
% m.walker@soton.ac.uk
% i.kuprov@soton.ac.uk

function hcch_tocsy_gb1()

% Protein data import
options.pdb_mol=1;
options.noshift='delete';
options.select='all';
[sys,inter]=protein('2N9K.pdb','2N9K.bmrb',options);

% Magnet field
sys.magnet=14.1;

% Algorithmic options
sys.tols.inter_cutoff=20.0;
sys.tols.prox_cutoff=4.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.level=4; bas.space_level=1;

% Algorithmic options
sys.enable={'greedy','prop_cache'};

% Sequence parameters
parameters.J_ch=140;
parameters.delta=1.1e-3;
parameters.sl_tmix=2e-3;
parameters.lamp=1e4;
parameters.dipsi_dur=22.5e-3;
parameters.sweep=[6000 13000 6000];
parameters.spins={'1H','13C','1H'};
parameters.offset=[2500 7000 2500];
parameters.npoints=[128 128 128];
parameters.zerofill=[256 256 256];
parameters.decouple_f3={'13C'};
parameters.axis_units='ppm';

% Create the spin system structure
spin_system=create(sys,inter);

% Kill nitrogens (not relevant)
spin_system=kill_spin(spin_system,strcmp('15N',spin_system.comp.isotopes));

% Build the basis
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@hcch_tocsy,parameters,'nmr');

% Apodization
fid=apodization(fid,'cosbell-3d');

% Fourier transform
spectrum=fftshift(fftn(fid,parameters.zerofill));

% Plotting
figure(); plot_3d(spin_system,abs(spectrum),parameters,...
                  10,[0.05 0.25 0.05 0.25],2,'positive');

end

