% 13C MAS spectrum of alanine powder (assuming decoupling of 1H),
% computed using the Floquet MAS formalism.
%
% Calculation time: minutes
%
% i.kuprov@soton.ac.uk

function mas_powder_ala_floquet()

% Spin system properties (PCM DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/alanine.log'),...
                    {{'C','13C'},{'N','15N'}},[182.1 264.5],[]);
% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.longitudinals={'15N'};
bas.projections=+1;

% Algorithmic options
sys.tols.inter_cutoff=5.0;
sys.tols.prox_cutoff=4.0;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.rate=2000;
parameters.axis=[1 1 1];
parameters.max_rank=17;
parameters.sweep=5e4;
parameters.npoints=256;
parameters.zerofill=1024;
parameters.offset=15000;
parameters.spins={'13C'};
parameters.decouple={};
parameters.axis_units='ppm';
parameters.invert_axis=1;
parameters.grid='rep_2ang_100pts_sph';
parameters.rho0=state(spin_system,'L+','13C','cheap');
parameters.coil=state(spin_system,'L+','13C','cheap');
parameters.verbose=0;

% Simulation
fid=floquet(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

