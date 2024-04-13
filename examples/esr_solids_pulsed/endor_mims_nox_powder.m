% Mims ENDOR simulation for a nitroxide radical powder. Ideal
% hard pulses are assumed.
%
% Calculation time: minutes
%
% i.kuprov@soton.ac.uk

function endor_mims_nox_powder()

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
parameters.offset=0;
parameters.npoints=512;
parameters.sweep=10e8;
parameters.tau=100e-9;
parameters.zerofill=4096;
parameters.spins={'E'};
parameters.axis_units='MHz';
parameters.grid='rep_2ang_12800pts_sph';

% Simulation
fid=powder(spin_system,@endor_mims,parameters,'esr');

% Crude apodization
fid=apodization(fid-mean(fid),'exp-1d',6);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);
kxlabel('Nuclear frequency, MHz');

end

