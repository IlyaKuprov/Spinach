% A deliberately adversarial example from Art Bochevarov at
% Schodinger Inc. In this case, IK-2 approximation in Liou-
% ville space generates an exceedingly large basis set; the
% calculation must instead be performed in Hilbert space 
% with permutation symmetry factorisation.
%
% Calculation time: hours.
%
% art.bochevarov@schrodinger.com
% ilya.kuprov@weizmann.ac.il

function fluoroisooctane()

% Magnet induction 
sys.magnet=11.74;

% Isotopes
sys.isotopes={'1H','1H','1H','1H','1H','1H',...
              '1H','1H','1H','1H','19F','1H',...
              '1H','1H','1H','1H','1H','1H'};

% Chemical shifts
inter.zeeman.scalar={0.994, 0.994, 0.994, ...
                     0.994, 0.994, 0.994, ...
                     0.994, 0.994, 0.994, ...
                     4.165, 0.000, 1.938, ...
                     1.039, 1.039, 1.039, ...
                     1.062, 1.062, 1.062};

% Larger J-couplings
inter.coupling.scalar=cell(18,18);
inter.coupling.scalar{10,11}=48.2;
inter.coupling.scalar{11,12}=23.6;
inter.coupling.scalar{10,12}= 6.5;
inter.coupling.scalar{12,13}= 6.0;
inter.coupling.scalar{12,14}= 6.0;
inter.coupling.scalar{12,15}= 6.0;
inter.coupling.scalar{12,16}= 6.0;
inter.coupling.scalar{12,17}= 6.0;
inter.coupling.scalar{12,18}= 6.0;

% Smaller J-couplings
inter.coupling.scalar{1,11}= 1.0;
inter.coupling.scalar{2,11}= 1.0;
inter.coupling.scalar{3,11}= 1.0;
inter.coupling.scalar{4,11}= 1.0;
inter.coupling.scalar{5,11}= 1.0;
inter.coupling.scalar{6,11}= 1.0;
inter.coupling.scalar{7,11}= 1.0;
inter.coupling.scalar{8,11}= 1.0;
inter.coupling.scalar{9,11}= 1.0;

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';
bas.sym_group={'S3','S3','S3'};
bas.sym_spins={[1 2 3],[4 5 6],[7 8 9]};

% Spinach housekeeping 
spin_system=create(sys,inter); 
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=1290;
parameters.sweep=2000;
parameters.npoints=4096;
parameters.zerofill=8192;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'gauss',6}});

% Fourier transform
spectrum=real(fftshift(fft(fid,parameters.zerofill)));

% Plotting
plot_1d(spin_system,spectrum,parameters);

end

