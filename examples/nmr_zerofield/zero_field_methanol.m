% Zero-field NMR spectroscopy - 13C methanol. Set to reproduce
% Figure 1 from http://dx.doi.org/10.1016/j.cplett.2013.06.042
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function zero_field_methanol()

% Spin system
sys.isotopes={'1H','1H','1H','13C'};

% Interactions
sys.magnet=0;
inter.zeeman.scalar={0 0 0 0};
inter.coupling.scalar{1,4}=141.0;
inter.coupling.scalar{2,4}=141.0;
inter.coupling.scalar{3,4}=141.0;
inter.coupling.scalar{4,4}=0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.sym_group={'S3'};
bas.sym_spins={[1 2 3]};

% Sequence parameters
parameters.sweep=700;
parameters.npoints=4096;
parameters.zerofill=16384;
parameters.offset=0;
parameters.spins={'1H'};
parameters.axis_units='Hz';
parameters.invert_axis=0;
parameters.flip_angle=pi/2;
parameters.detection='uniaxial';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@zerofield,parameters,'labframe');

% Apodisation
fid=apodisation(spin_system,fid-mean(fid),{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters); 

end

