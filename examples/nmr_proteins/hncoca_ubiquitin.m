% Theoretical HN(CO)CA of human ubiquitin. It is assumed that 
% only the backbone is 13C,15N-labelled.
%
% Calculation time: minutes, faster with a Tesla A100 GPU.
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk

function hncoca_ubiquitin()

% Protein data import
options.pdb_mol=1;
options.noshift='delete';
options.select='backbone-minimal';
[sys,inter]=protein('1D3Z.pdb','1D3Z.bmrb',options);

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
parameters.spins={'15N','13C','1H'};
parameters.offset=[-7100 8450 4850];
parameters.sweep=[2500 4500 3000];
parameters.npoints=[64 64 64];
parameters.zerofill=[256 256 256];
parameters.axis_units='ppm';

% Simulation
fid=liquid(spin_system,@hncoca,parameters,'nmr');

% Apodization
fid=apodization(fid,'sqcosbell-3d');

% Fourier transform
spectrum=fftn(fid,parameters.zerofill);
spectrum=fftshift(spectrum);

% Plotting
figure(); plot_3d(spin_system,abs(spectrum),parameters,...
                  10,[0.2 0.9 0.2 0.9],2,'positive');

end

