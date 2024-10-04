% Zero-field NMR spectroscopy - 15N pyridine. Set to reproduce
% Fig 2 from http://dx.doi.org/10.1103/PhysRevLett.107.107601
%
% Calculation time: seconds
%
% ortmeier@ncsu.edu

function small_field_formic_acid()

% Magnetic field
sys.magnet=1.76e-7;

% Spin system
sys.isotopes={'1H','13C'};

% Interactions
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=221; 

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Sequence parameters
parameters.sweep=600;
parameters.npoints=8196;
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
fid=apodisation(spin_system,fid-mean(fid),{{'exp',12}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

