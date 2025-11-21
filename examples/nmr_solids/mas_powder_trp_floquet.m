% 13C MAS spectrum of tryptophan powder (assuming decoupling of 1H),
% computed using the Floquet MAS formalism. Isotropic chemical shifts
% come from the experimental data. Coordinates are from X-ray data 
% and CSAs are estimated with DFT.
%
% Calculation time: days (hours with a Tesla card)
%
% ilya.kuprov@weizmann.ac.il
% giuseppe.pileio@soton.ac.uk
% maria.concistre@soton.ac.uk

function mas_powder_trp_floquet()

%% First molecule in the unit cell

% Spin system properties (DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/trp_xray.out'),...
                         {{'C','13C'},{'N','15N'}},[174.0 0],[]);
% Magnet field
sys.magnet=9.4;

% First conformation
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

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.rate=14000;
parameters.axis=[1 1 1];
parameters.max_rank=11;
parameters.grid='leb_2ang_rank_11';
parameters.sweep=1e5;
parameters.npoints=2048;
parameters.zerofill=8192;
parameters.offset=0;
parameters.spins={'13C'};
parameters.decouple={};
parameters.axis_units='ppm';
parameters.invert_axis=1;
parameters.rho0=state(spin_system,'L+','13C','cheap');
parameters.coil=state(spin_system,'L+','13C','cheap');
parameters.verbose=1;

% Simulation
fid=floquet(spin_system,@acquire,parameters,'nmr');

%% Second molecule in the unit cell

% Spin system properties (DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/trp_xray.out'),...
                         {{'C','13C'},{'N','15N'}},[174.0 0],[]);
% Magnet field
sys.magnet=9.4;

% Algorithmic options
sys.tols.inter_cutoff=5.0;
sys.tols.prox_cutoff=4.0;
sys.disable={'trajlevel'};

% Second conformation
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,2,124.2);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,3,110.1);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,6,114.7);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,4,118.0);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,5,119.3);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,7,107.5);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,8,134.9);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,9,125.0);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,10,28.0);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,11,52.1);
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,12,173.3);

% Simulation
fid=fid+floquet(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

