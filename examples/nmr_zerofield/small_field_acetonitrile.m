% Small-field NMR spectroscopy - acetonitrile with 13C on the 
% methyl group. Set to reproduce Figure 3 from 
%
%       http://dx.doi.org/10.1103/PhysRevLett.107.107601
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function small_field_acetonitrile()

% Magnetic field, 2.64 mG
sys.magnet=2.64e-3*1e-4;

% Spin system
sys.isotopes={'1H','1H','1H','13C'};

% Interactions
inter.coupling.scalar=cell(4,4);
inter.coupling.scalar{1,4}=136.200;
inter.coupling.scalar{2,4}=136.200;
inter.coupling.scalar{3,4}=136.200;

% Temperature
inter.temperature=298;

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

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

