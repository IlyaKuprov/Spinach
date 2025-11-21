% A soft pulse simulation for a nitroxide radical powder. The soft 
% pulse is simulated using the Fokker-Planck formalism; it is fol-
% lowed by time domain acquisition and Fourier transform.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function spa_nitroxide_powder()

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
parameters.pulse_phi=-pi/2;
parameters.pulse_frq=-300e6;
parameters.pulse_dur=100e-9;
parameters.pulse_pwr=2*pi*16.5e6;
parameters.method='expm';

% Simulation
fid=powder(spin_system,@sp_acquire,parameters,'esr');

% Apodisation
fid=apodisation(spin_system,fid,{{'crisp'}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

