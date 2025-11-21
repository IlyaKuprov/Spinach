% Zero-field NMR spectroscopy - acetonitrile. The simulation 
% proceeds by computing the exact thermal equilibrium state 
% and them propagating it through a time-dependent field drop.
% Set to reproduce Figure 7 from 
%
%          https://doi.org/10.1016/j.jmr.2017.08.016
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function field_drop_acetonitrile()

% Magnetic field (polariser)
sys.magnet=2.0;

% Spin system
sys.isotopes={'1H','1H','1H','13C','13C','15N'};

% Interactions
inter.coupling.scalar{1,4}=136.200;
inter.coupling.scalar{2,4}=136.200;
inter.coupling.scalar{3,4}=136.200;
inter.coupling.scalar{1,5}=-9.924;
inter.coupling.scalar{2,5}=-9.924;
inter.coupling.scalar{3,5}=-9.924;
inter.coupling.scalar{1,6}=-1.688;
inter.coupling.scalar{2,6}=-1.688;
inter.coupling.scalar{3,6}=-1.688;
inter.coupling.scalar{4,5}=57.010;
inter.coupling.scalar{4,6}=2.822;
inter.coupling.scalar{5,6}=-17.419;
inter.coupling.scalar{6,6}=0;

% Temperature
inter.temperature=298;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.sym_group={'S3'};
bas.sym_spins={[1 2 3]};

% Sequence parameters
parameters.sweep=700;
parameters.npoints=4096;
parameters.zerofill=16586;
parameters.offset=0;
parameters.spins={'1H'};
parameters.axis_units='Hz';
parameters.invert_axis=0;
parameters.drop_field=0.1e-3;  % Tesla
parameters.drop_rate=10;       % Hz
parameters.drop_time=0.5;      % seconds
parameters.drop_npoints=100;
parameters.flip_angle=pi/2;
parameters.detection='uniaxial';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@zulf_abrupt,parameters,'labframe');

% Apodisation
fid=apodisation(spin_system,fid-mean(fid),{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

