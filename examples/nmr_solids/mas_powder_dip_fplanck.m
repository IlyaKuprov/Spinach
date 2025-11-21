% Spinning powder pulse-acquire experiment on a two-spin system 
% with a dipolar coupling using Fokker-Planck formalism:
%
%           https://doi.org/10.1016/j.jmr.2016.07.005
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function mas_powder_dip_fplanck()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={5.0,-2.0};
inter.coordinates={[0.0 0.0 0.0];
                   [0.0 3.9 0.1]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=+1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Pulse-acquire setup
parameters.rate=1000;
parameters.axis=[1 1 1];
parameters.max_rank=17;
parameters.grid='leb_2ang_rank_17';
parameters.sweep=2e4;
parameters.npoints=512;
parameters.zerofill=4096;
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H','cheap');
parameters.coil=state(spin_system,'L+','1H','cheap');
parameters.verbose=0;

% Simulation
fid=singlerot(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

