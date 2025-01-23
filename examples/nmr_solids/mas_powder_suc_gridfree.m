% 13C MAS spectrum of sucrose powder (assuming decoupling of 1H),
% computed using the grid-free Fokker-Planck MAS formalism. Che-
% mical shielding tensors, J-couplings and coordinates are esti-
% mated with DFT. A polyadic representation of the evolution ge-
% nerator is used, further particulars here:
%
%             https://doi.org/10.1126/sciadv.aaw8962
%
% Calculation time: hours on a Tesla V100 GPU,
%                   much longer on CPU
%
% ilya.kuprov@weizmann.ac.il

function mas_powder_suc_gridfree()

% Spin system properties (PCM DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/sucrose.log'),...
                                        {{'C','13C'}},182.1,[]);        
% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.projections=+1;
bas.level=3;

% Algorithmic options
sys.tols.inter_cutoff=5.0;
sys.tols.prox_cutoff=4.0;
sys.enable={'greedy','polyadic','gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.axis=[1 1 1];
parameters.max_rank=23;
parameters.rate=6000;
parameters.sweep=5e4;
parameters.npoints=256;
parameters.zerofill=1024;
parameters.offset=15000;
parameters.spins={'13C'};
parameters.decouple={};
parameters.axis_units='ppm';
parameters.invert_axis=1;
parameters.rho0=state(spin_system,'L+','13C','cheap');
parameters.coil=state(spin_system,'L+','13C','cheap');
parameters.verbose=1;

% Run the simulation
fid=gridfree(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

