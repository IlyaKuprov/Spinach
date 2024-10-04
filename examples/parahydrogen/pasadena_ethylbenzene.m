% PASADENA experiment simulation for the parahydrogenation of styrene
% into ethylbenzene. Set to reproduce the top trace of Fig 5 in
%
%                https://doi.org/10.1039/b914188j
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function pasadena_ethylbenzene()

% Spin system
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'};

% Magnetic field
sys.magnet=7.05;

% Chemical shifts
inter.zeeman.scalar={1.201 1.201 1.201 2.625 2.625...
                     7.207 7.207 7.265 7.265 7.155};

% Scalar couplings
inter.coupling.scalar=cell(10);
inter.coupling.scalar{1,4}=7.63;   inter.coupling.scalar{2,4}=7.63;
inter.coupling.scalar{3,4}=7.63;   inter.coupling.scalar{1,5}=7.63;
inter.coupling.scalar{2,5}=7.63;   inter.coupling.scalar{3,5}=7.63;
inter.coupling.scalar{6,8}=7.63;   inter.coupling.scalar{6,10}=1.26;
inter.coupling.scalar{7,9}=7.63;   inter.coupling.scalar{7,10}=1.26;
inter.coupling.scalar{8,10}=7.44;  inter.coupling.scalar{9,10}=7.44;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;
bas.sym_group={'S3','S2'};
bas.sym_spins={[1 2 3],[4 5]};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,{'Lz','Lz'},{1,4});
parameters.coil=state(spin_system,'L+','1H');
parameters.pulse_op=(operator(spin_system,'L+','1H')-...
                     operator(spin_system,'L-','1H'))/2i;
parameters.pulse_angle=-pi/4;
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
fid=apodisation(spin_system,fid,{{'gauss',10}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

