% A hole burning simulation for a nitroxide radical. The soft pulse 
% is simulated using Fokker-Planck formalism; it is followed by an
% ideal hard pulse, acquisition and Fourier transform.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function holeburn_nitroxide_powder()

% Isotopes                          
sys.isotopes={'E','14N'};
                          
% Magnet field
sys.magnet=3.5;

% Interactions
inter.zeeman.matrix=cell(1,2);
inter.zeeman.matrix{1}=[2.01045  0.00000  0.00000
                        0.00000  2.00641  0.00000
                        0.00000  0.00000  2.00211];
inter.coupling.matrix=cell(2,2);
inter.coupling.matrix{1,2}=[1.2356  0.0000  0.6322
                            0.0000  1.1266  0.0000
                            0.6322  0.0000  8.2230]*1e7;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Disable trajectory-level SSR algorithms
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'Lz','E');
parameters.coil=state(spin_system,'L+','E');
parameters.decouple={};
parameters.offset=-2e8;
parameters.sweep=8e8;
parameters.npoints=64;
parameters.zerofill=512;
parameters.axis_units='MHz';
parameters.grid='rep_2ang_3200pts_sph';
parameters.derivative=0;
parameters.invert_axis=0;

% Soft pulse parameters
parameters.pulse_rnk=2;
parameters.pulse_phi=pi/2;
parameters.method='expv';
parameters.pulse_dur=1e-9*ones(100,1);

% Simulation A - chirp
parameters.pulse_pwr=2*pi*10e6*ones(100,1);
parameters.pulse_frq=linspace(-350e6,-250e6,100);
fid_a=powder(spin_system,@holeburn,parameters,'esr');

% Simulation B - fixed frequency
parameters.pulse_pwr=2*pi*10e6*ones(100,1);
parameters.pulse_frq=-300e6*ones(100,1);
fid_b=powder(spin_system,@holeburn,parameters,'esr');

% Simulation C - reference
parameters.pulse_pwr=zeros(100,1);
parameters.pulse_frq=zeros(100,1);
fid_c=powder(spin_system,@holeburn,parameters,'esr');

% Apodization
fid_a=apodization(fid_a,'exp-1d',6);
fid_b=apodization(fid_b,'exp-1d',6);
fid_c=apodization(fid_c,'exp-1d',6);

% Fourier transform
spectrum_a=fftshift(fft(fid_a,parameters.zerofill));
spectrum_b=fftshift(fft(fid_b,parameters.zerofill));
spectrum_c=fftshift(fft(fid_c,parameters.zerofill));

% Plotting
figure(); hold on;
plot_1d(spin_system,real(spectrum_a),parameters,'r-'); 
plot_1d(spin_system,real(spectrum_b),parameters,'b-');
plot_1d(spin_system,real(spectrum_c),parameters,'k-');
legend({'chirp pulse','soft pulse','reference'});

end

