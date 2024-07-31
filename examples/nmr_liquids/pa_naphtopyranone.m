% NMR spectrum of 3-phenylmethylene-1H,3H-naphtho-[1,8-c,d]-pyran-1-one,
% magnetic parameters from:
%
%              http://dx.doi.org/10.1016/j.saa.2010.11.015
%
% Calculation time: seconds
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk

function pa_naphtopyranone()

% Magnetic induction
sys.magnet=14.095;

% Spin system
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'};

% Chemical shifts
inter.zeeman.scalar={8.345,7.741,8.097,8.354,7.784,8.330,7.059,...
                     7.941,7.466,7.326,7.466,7.941};

% Scalar couplings
inter.coupling.scalar=cell(12,12);
inter.coupling.scalar{1,2}=7.8;
inter.coupling.scalar{1,3}=0.9;
inter.coupling.scalar{2,3}=7.8;
inter.coupling.scalar{4,5}=8.4;
inter.coupling.scalar{4,6}=1.2;
inter.coupling.scalar{5,6}=7.2;
inter.coupling.scalar{8,9}=7.8;
inter.coupling.scalar{8,10}=1.2;
inter.coupling.scalar{9,10}=7.8;
inter.coupling.scalar{10,11}=7.8;
inter.coupling.scalar{10,12}=1.2;
inter.coupling.scalar{11,12}=7.8;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=4600;
parameters.sweep=1200;
parameters.npoints=4096;
parameters.zerofill=32768;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodization
fid=apodization(fid,'exp-1d',10);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

