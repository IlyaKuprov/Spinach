% Zero-field NMR spectroscopy - 15N pyridine. Set to reproduce
% Figure 3 from http://dx.doi.org/10.1021/ja2112405
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function zero_field_pyridine()

% Magnetic field
sys.magnet=0;

% Spin system
sys.isotopes={'1H','1H','1H','1H','1H','15N'};

% Interactions
inter.coupling.scalar{1,2}=   4.88; 
inter.coupling.scalar{4,5}=   4.88;
inter.coupling.scalar{1,4}=   0.97; 
inter.coupling.scalar{2,5}=   0.97;
inter.coupling.scalar{1,3}=   1.83; 
inter.coupling.scalar{3,5}=   1.83;
inter.coupling.scalar{1,5}=  -0.12;
inter.coupling.scalar{2,3}=   7.62; 
inter.coupling.scalar{3,4}=   7.62;
inter.coupling.scalar{2,4}=   1.38;
inter.coupling.scalar{1,6}= -10.93; 
inter.coupling.scalar{5,6}= -10.93;
inter.coupling.scalar{2,6}=  -1.47; 
inter.coupling.scalar{4,6}=  -1.47;
inter.coupling.scalar{3,6}=   0.27;
inter.coupling.scalar{6,6}=   0.00;

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Sequence parameters
parameters.sweep=60;
parameters.npoints=512;
parameters.zerofill=1024;
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
fid=apodisation(spin_system,fid-mean(fid),{{'exp',12}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters); 

end

