% Constant-time HSQC experiment simulation for the 
% GB1 protein.
%
% Simulation time: hours, faster with a Tesla A100 GPU.
%
% m.walker@soton.ac.uk
% i.kuprov@soton.ac.uk

function ct_hsqc_gb1()

% Protein data import
options.pdb_mol=1;
options.noshift='delete';
options.select='backbone';
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

% Sequence parameters
parameters.J=90;
parameters.sweep=[3000 3000];
parameters.offset=[-7300 5100];
parameters.npoints=[128 128];
parameters.zerofill=[512 512];
parameters.spins={'15N','1H'};
parameters.decouple_f2={'15N'};
parameters.axis_units='ppm';

% Create the spin system structure
spin_system=create(sys,inter);

% Kill nitrogens (protein assumed unlabelled)
spin_system=kill_spin(spin_system,strcmp('13C',spin_system.comp.isotopes));

% Build the basis
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@ct_hsqc,parameters,'nmr');

% Apodization
fid.pos=apodization(fid.pos,'sqcosbell-2d');
fid.neg=apodization(fid.neg,'sqcosbell-2d');

% F2 Fourier transform
f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1);
f1_neg=fftshift(fft(fid.neg,parameters.zerofill(2),1),1);

% Form States signal
fid=f1_pos+conj(f1_neg);

% F1 Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill(1),2),2);

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,abs(spectrum),parameters,...
        20,[0.1 0.5 0.1 0.5],2,256,6,'positive');

end

