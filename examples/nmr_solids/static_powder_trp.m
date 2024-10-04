% 13C NMR spectrum of tryptophan powder. Isotropic chemical shifts
% come from the experimental data. Coordinates and CSAs are estima-
% ted with DFT. Protons are assumed to be decoupled.
%
% Calculation time: hours
%
% i.kuprov@soton.ac.uk

function static_powder_trp()

% Spin system properties (DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/trp_xray.out'),...
                         {{'C','13C'},{'N','15N'}},[174.0 0],[]);
% Magnet field
sys.magnet=14.1;

% Experimental chemical shifts
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,2,124.2);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,3,110.1);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,6,114.7);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,4,118.0);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,5,119.3);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,7,107.5);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,8,134.9);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,9,125.0);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,10,26.8);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,11,54.6);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,12,174.4);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.longitudinals={'15N'};
bas.projections=+1;
bas.level=3;

% Algorithmic options
sys.tols.inter_cutoff=5.0;
sys.tols.prox_cutoff=4.0;
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.sweep=6e4;
parameters.npoints=128;
parameters.zerofill=512;
parameters.offset=18000;
parameters.spins={'13C'};
parameters.decouple={};
parameters.axis_units='ppm';
parameters.invert_axis=1;
parameters.grid='rep_2ang_6400pts_sph';
parameters.rho0=state(spin_system,'L+','13C','cheap');
parameters.coil=state(spin_system,'L+','13C','cheap');
parameters.verbose=0;

% Simulation
fid=powder(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

