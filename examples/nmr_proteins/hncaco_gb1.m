% Simulated HNCACO spectrum of GB1 protein. It is assumed that
% only the backbone is 13C,15N-labelled.
%
% Calculation time: minutes, faster with a Tesla A100 GPU.
%
% m.walker@soton.ac.uk
% i.kuprov@soton.ac.uk

function hncaco_gb1()

% Protein data import
options.pdb_mol=1;
options.noshift='delete';
options.select='backbone-minimal';
[sys,inter]=protein('2N9K.pdb','2N9K.bmrb',options);

% Magnet field
sys.magnet=14.1;

% Tolerances
sys.tols.inter_cutoff=2.0;
sys.tols.prox_cutoff=4.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.level=4; bas.space_level=1;

% Algorithmic options
sys.enable={'greedy','gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.J_nh=92;
parameters.T=25e-3;
parameters.delta2=3e-3;
parameters.spins={'15N','13C','1H'};
parameters.sweep=[3000 2500 3000];
parameters.offset=[-7200 26500 5100];
parameters.npoints=[64 64 64];
parameters.zerofill=[256 256 256];
parameters.axis_units='ppm';

% Simulation
fid=liquid(spin_system,@hncaco,parameters,'nmr');

% Apodization
fid=apodization(fid,'cosbell-3d');

% Fourier transform
spectrum=fftshift(fftn(fid,parameters.zerofill));

% Plotting
figure(); plot_3d(spin_system,abs(spectrum),parameters,...
                  10,[0.2 0.9 0.2 0.9],2,'positive');

end

