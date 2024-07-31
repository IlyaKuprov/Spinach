% Constant-time COSY experiment simulation for the 
% GB1 protein.
%
% Simulation time: minutes, faster with a Tesla A100 GPU.
%
% m.walker@soton.ac.uk
% i.kuprov@soton.ac.uk

function ct_cosy_gb1()

% Protein data import
options.pdb_mol=1;
options.noshift='delete';
options.select='backbone';
[sys,inter]=protein('2N9K.pdb','2N9K.bmrb',options);

% Magnet field
sys.magnet=14.1;

% Tolerances
sys.tols.inter_cutoff=1.0;
sys.tols.prox_cutoff=3.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.level=4; bas.space_level=1;

% Algorithmic options
sys.enable={'greedy','gpu'};

% Sequence parameters
parameters.offset=2400;
parameters.sweep=[9000 9000];
parameters.npoints=[256 256];
parameters.zerofill=[512 512];
parameters.spins={'1H'};
parameters.axis_units='ppm';

% Spinach housekeeping
spin_system=create(sys,inter);

% Kill carbons and nitrogens (protein assumed unlabelled)
spin_system=kill_spin(spin_system,strcmp('13C',spin_system.comp.isotopes));
spin_system=kill_spin(spin_system,strcmp('15N',spin_system.comp.isotopes));

% Build the basis
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@ct_cosy,parameters,'nmr');

% Apodization
fid=apodization(fid,'cosbell-2d');

% Fourier transform
spectrum=fftn(fid,parameters.zerofill);
spectrum=fftshift(spectrum);

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,abs(spectrum),parameters,...
        20,[0.01 0.125 0.01 0.125],2,256,6,'positive');

end

