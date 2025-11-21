% ALTADENA experiment simulation for the parahydrogenation of acrolein
% into propanal. Simple model of the ALTADENA effect is used: perfect-
% ly adiabatic transfer is assumed and the isotropic mixing in low fi-
% eld is ignored completely. Note the small flip angle.
%
% Calculation time: seconds
%
% Ronghui Zhou (hui@ufl.edu)
% ilya.kuprov@weizmann.ac.il

function altadena_propanal()

% Spin system
sys.isotopes={'1H','1H','1H','1H','1H','1H'};

% Magnetic field
sys.magnet=7.05;

% Chemical shifts
inter.zeeman.scalar={1.11 1.11 1.11 2.46 2.46 9.79};

% Scalar couplings
inter.coupling.scalar=cell(6);
inter.coupling.scalar{1,4}=7.3;   inter.coupling.scalar{2,4}=7.3;
inter.coupling.scalar{3,4}=7.3;   inter.coupling.scalar{1,5}=7.3;
inter.coupling.scalar{2,5}=7.3;   inter.coupling.scalar{3,5}=7.3;
inter.coupling.scalar{4,6}=1.4;   inter.coupling.scalar{5,6}=1.4;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.sym_group={'S3','S2'};
bas.sym_spins={[1 2 3],[4 5]};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,{'Lz','Lz'},{1,4})-...
                0.5*state(spin_system,{'Lz'},{1})+...
                0.5*state(spin_system,{'Lz'},{4});
parameters.coil=state(spin_system,'L+','1H');
parameters.pulse_op=(operator(spin_system,'L+','1H')-...
                     operator(spin_system,'L-','1H'))/2i;
parameters.pulse_angle=pi/100;
parameters.decouple={};
parameters.offset=500;
parameters.sweep=1000;
parameters.npoints=1024;
parameters.zerofill=8192;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@hp_acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

