% Zero-field NMR spectroscopy - benzene with one 13C 
% nucleus. Set to reproduce Figure 2 from 
%     
%     https://doi.org/10.1016/j.jmr.2017.08.016
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function zero_field_benzene()

% Magnetic field
sys.magnet=0;

% Spin system
sys.isotopes={'1H','1H','1H','1H','1H','1H','13C'};

% Interactions
inter.coupling.scalar{1,7}=158.363;
inter.coupling.scalar{2,7}=1.136;
inter.coupling.scalar{3,7}=7.609;
inter.coupling.scalar{4,7}=-1.285;
inter.coupling.scalar{5,7}=7.609;
inter.coupling.scalar{6,7}=1.136;
inter.coupling.scalar{1,2}=7.534;
inter.coupling.scalar{1,3}=1.381;
inter.coupling.scalar{1,4}=0.658;
inter.coupling.scalar{1,5}=1.381;
inter.coupling.scalar{1,6}=7.534;
inter.coupling.scalar{2,3}=7.543;
inter.coupling.scalar{2,4}=1.382;
inter.coupling.scalar{2,5}=0.660;
inter.coupling.scalar{2,6}=1.384;
inter.coupling.scalar{3,4}=7.543;
inter.coupling.scalar{3,5}=1.387;
inter.coupling.scalar{3,6}=0.660;
inter.coupling.scalar{4,5}=7.543;
inter.coupling.scalar{4,6}=1.382;
inter.coupling.scalar{5,6}=7.543;
inter.coupling.scalar{7,7}=0;

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Sequence parameters
parameters.sweep=400;
parameters.npoints=8162;
parameters.zerofill=16586;
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

