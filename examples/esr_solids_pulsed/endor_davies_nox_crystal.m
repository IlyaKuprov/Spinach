% Davies ENDOR simulation for a nitroxide radical at a single
% orientation. Soft pulses are simulated using Fokker-Planck 
% formalism. 
%
% Calculation time: minutes
%
% i.kuprov@soton.ac.uk

function endor_davies_nox_crystal()

% Isotopes                          
sys.isotopes={'E','14N'};
                          
% Magnet field
sys.magnet=3.35;

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

% Relaxation theory
inter.relaxation={'t1_t2'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.r1_rates={20e6 0.5e6};
inter.r2_rates={20e6 0.5e6};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

%% Stage 1: pulse-acquire ESR spectrum

% Sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'L+','E');
parameters.coil=state(spin_system,'L+','E');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=5e8;
parameters.npoints=128;
parameters.zerofill=512;
parameters.axis_units='MHz';
parameters.orientation=[0 0 0];
parameters.derivative=0;
parameters.invert_axis=1;

% Simulation
fid=crystal(spin_system,@acquire,parameters,'esr');

% Apodisation
fid=apodisation(spin_system,fid,{{'none'}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); scale_figure([1.5 1.5]); subplot(2,2,1);
plot_1d(spin_system,real(spectrum),parameters);
title('Single crystal ESR spectrum');

%% Stage 2: ENDOR, left line

% Sequence parameters
clear('parameters');
parameters.spins={'E','14N'};
parameters.rho0=state(spin_system,'Lz','E');
parameters.coil=state(spin_system,'L+','E');
parameters.orientation=[0 0 0];
parameters.offset=[0 0];
parameters.method='expm';

% Electron pulse parameters
parameters.e_rnk=2;
parameters.e_phi=0;
parameters.e_dur=10e-9;
parameters.e_pwr=2*pi*16.5e7;

% Nucleus pulse parameters
parameters.n_rnk=3;
parameters.n_phi=0;
parameters.n_frq=linspace(-200e6,200e6,200);
parameters.n_dur=1e-7;
parameters.n_pwr=pi*1e7;
 
% ENDOR simulation - three electron lines
parameters.e_frq=+93e6;
answer_a=crystal(spin_system,@endor_davies,parameters,'esr');
parameters.e_frq=+10e6;
answer_b=crystal(spin_system,@endor_davies,parameters,'esr');
parameters.e_frq=-73e6;
answer_c=crystal(spin_system,@endor_davies,parameters,'esr');

% Plotting
subplot(2,2,2); plot(parameters.n_frq/1e6,real(answer_a));
kgrid; axis tight; xlabel('Nuclear frequency, MHz');
ylabel('(RF on)/(RF off)'); title('Single crystal ENDOR');
subplot(2,2,3); plot(parameters.n_frq/1e6,real(answer_b));
kgrid; axis tight; xlabel('Nuclear frequency, MHz');
ylabel('(RF on)/(RF off)'); title('Single crystal ENDOR');
subplot(2,2,4); plot(parameters.n_frq/1e6,real(answer_c));
kgrid; axis tight; xlabel('Nuclear frequency, MHz');
ylabel('(RF on)/(RF off)'); title('Single crystal ENDOR');

end

