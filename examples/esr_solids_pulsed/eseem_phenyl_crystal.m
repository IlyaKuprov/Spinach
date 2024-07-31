% Two-pulse X-band ESEEM spectrum of a phenyl radical at a specific
% orientation relative to the lab frame. Magnetic parameters taken 
% from a DFT calculation. Ideal pulses are assumed.
%
% Calculation time: minutes
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk

function eseem_phenyl_crystal()

% Spin system properties (vacuum DFT calculation)
options.no_xyz=1;
[sys,inter]=g2spinach(gparse('../standard_systems/phenyl.log'),...
                      {{'E','E'},{'H','1H'}},[0 0],options);
% Magnet field
sys.magnet=0.33;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'Lz','E');
parameters.coil=state(spin_system,'L+','E');
parameters.screen=state(spin_system,'L-','E');
parameters.pulse_op=(operator(spin_system,'L+','E')-...
                     operator(spin_system,'L-','E'))/2i;
parameters.npoints=512;
parameters.timestep=1e-8;
parameters.orientation=[pi/5 pi/4 pi/3];
parameters.zerofill=4096;

% Simulation
fid=crystal(spin_system,@eseem,parameters,'esr');

% Plot the time domain signal
figure(); subplot(2,1,1);
plot((0:(parameters.npoints-1))*parameters.timestep*1e6,real(fid));
xlabel('time, \mus'); axis tight; kgrid;

% Crude apodization
fid=apodization(fid-mean(fid),'kaiser-1d',6);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));
ax=linspace(-1/(parameters.timestep),1/(parameters.timestep),...
            parameters.zerofill)*1e-6;

% Plot the spectrum
subplot(2,1,2); plot(ax,abs(spectrum)); 
xlabel('Frequency, MHz'); axis tight; kgrid;

end

