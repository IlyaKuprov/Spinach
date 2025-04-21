% 13C 2D PDSD spectrum of ubiquitin.
%
% Calculation time: hours.
%
% ilya.kuprov@weizmann.ac.il

function pdsd_ubiquitin()

% Protein data import
options.pdb_mol=1;
options.select=1:32; % 'all';
options.noshift='delete';
[sys,inter]=protein('1D3Z.pdb','1D3Z.bmrb',options);

% Magnet field
sys.magnet=21.1356;

% Tolerances
sys.tols.inter_cutoff=2.0;
sys.tols.prox_cutoff=5.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.level=3; bas.space_level=3;

% Algorithmic options
sys.enable={'greedy','gpu'};

% Create the spin system structure
spin_system=create(sys,inter);

% Remove all nitrogens atoms
n_idx=strcmp('15N',spin_system.comp.isotopes);
spin_system=kill_spin(spin_system,n_idx);

% Build the basis
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.spins={'1H','13C'};
parameters.rate=100000;
parameters.tmix=1e-3;
parameters.max_rank=3;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.grid='rep_2ang_100pts_sph';
parameters.npoints=[256 256];
parameters.zerofill=[1024 1024];
parameters.offset=[5000 20000];
parameters.sweep=50000;
parameters.verbose=1;

% Simulation
fid=singlerot(spin_system,@pdsd,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'cos'},{'cos'}});

% Fourier transform
spectrum=fftshift(fft2(fid,parameters.zerofill(2),...
                           parameters.zerofill(1)));

% Plotting
figure(); scale_figure([1.5 2.0]);
parameters.spins={'13C'};
parameters.offset=parameters.offset(2);
plot_2d(spin_system,abs(spectrum),parameters,...
        20,[0.01 0.1 0.01 0.1],2,256,6,'positive'); 

end

