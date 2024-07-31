% Powder-averaged two-pulse ESEEM on a 14N nitroxide radical. Time-domain
% simulation in Liouville space with powder averaging over a finite grid.
% Set to reproduce Figure 4a in http://dx.doi.org/10.1063/1.453532, ideal
% pulses are assumed.
%
% Calculation time: seconds
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk

function eseem_nitroxide_powder()

% Magnet field
sys.magnet=0.3249;

% System specification
sys.isotopes={'14N','E'};
inter.coupling.eigs=cell(2,2);
inter.coupling.euler=cell(2,2);
inter.coupling.eigs{1,1}=[-0.4 -1.6 +2.0]*1e5;
inter.coupling.eigs{1,2}=[2.0 2.0 2.0]*1e6;
inter.coupling.euler{1,1}=[0 0 0];
inter.coupling.euler{1,2}=[0 0 0];

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Disable trajectory-level SSR algorithms
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set the sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'Lz','E');
parameters.coil=state(spin_system,'L+','E');
parameters.screen=state(spin_system,'L-','E');
parameters.pulse_op=(operator(spin_system,'L+','E')-...
                     operator(spin_system,'L-','E'))/2i;
parameters.offset=0;
parameters.npoints=512;
parameters.timestep=2e-7;
parameters.zerofill=2048;
parameters.grid='rep_2ang_400pts_sph';

% Simulation
fid=powder(spin_system,@eseem,parameters,'esr');

% Run apodization
fid=apodization(mean(fid)-fid,'exp-1d',5);

% Run Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plot the time domain signal
figure(); subplot(2,1,1);
plot((0:(parameters.npoints-1))*parameters.timestep*1e6,real(fid));
xlabel('time, \mus'); axis tight; kgrid;

% Plot the spectrum
subplot(2,1,2);
plot(linspace(-1/parameters.timestep,1/parameters.timestep,...
     parameters.zerofill)*1e-6,real(spectrum));
xlabel('frequency, MHz'); axis tight; kgrid;

end

