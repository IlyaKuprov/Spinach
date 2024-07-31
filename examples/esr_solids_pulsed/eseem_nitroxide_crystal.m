% Two-pulse X-band ESEEM spectrum of a nitroxide radical at a specific
% orientation relative to the lab frame. Magnetic parameters taken from
% a DFT calculation. Ideal pulses are assumed.
%
% Calculation time: seconds
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk

function eseem_nitroxide_crystal()

% Isotopes                          
sys.isotopes={'E','14N'};
                          
% Interactions
inter.zeeman.matrix=cell(1,2);
inter.zeeman.matrix{1}=[2.01045  0.00000  0.00000
                        0.00000  2.00641  0.00000
                        0.00000  0.00000  2.00211];
inter.coupling.matrix=cell(2,2);
inter.coupling.matrix{1,2}=[1.2356  0.0000  0.6322
                            0.0000  1.1266  0.0000
                            0.6322  0.0000  8.2230]*1e7;

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
parameters.npoints=1024;
parameters.timestep=1.25e-8;
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

% Plot the spectrum
subplot(2,1,2);
plot(linspace(-1/(parameters.timestep),1/(parameters.timestep),...
     parameters.zerofill)*1e-6,abs(spectrum));
xlabel('frequency, MHz'); axis tight; kgrid;

end

