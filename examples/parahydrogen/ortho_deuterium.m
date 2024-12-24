% Ortho-deuteration simulation for acrylonitrile in Figure 1 of
% the paper by Natterer, Greve, and Bargon:
%
%         https://doi.org/10.1016/S0009-2614(98)00784-2
%
% Simulation time: seconds
%
% ilya.kuprov@weizmann.ac.il
% anakin.aden@mpinat.mpg.de 

function ortho_deuterium()

% Bargon's magnet
sys.magnet=4.697;

% Deuteration product
sys.isotopes={'1H','1H','2H','1H','2H'};
inter.zeeman.scalar={1.2 1.2 1.2 2.3 2.3};
inter.coupling.scalar=cell(5,5);
inter.coupling.scalar{1,3}=2.0;
inter.coupling.scalar{2,3}=2.0;
inter.coupling.scalar{4,5}=2.0;
inter.coupling.scalar{3,4}=1.2;
inter.coupling.scalar{1,5}=1.2;
inter.coupling.scalar{2,5}=1.2;
inter.coupling.scalar{3,5}=0.2;

% Hilbert space
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Continuous deuteration
options.dephasing=1;

% Singlet and quintet on deuterium 
[S,~,Q]=deut_pair(spin_system,3,5,options);

% Experiment parameters
parameters.spins={'2H'};
parameters.rho0=S+Q{1}+Q{2}+Q{3}+Q{4}+Q{5};
parameters.coil=state(spin_system,'L+','2H');
parameters.pulse_op=operator(spin_system,'Ly','2H');
parameters.pulse_angle=pi/4;
parameters.decouple={};
parameters.offset=50;
parameters.sweep=120;
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.axis_units='ppm';

% Simulation
fid=liquid(spin_system,@hp_acquire,parameters,'nmr');

% Apodisation and sign flip
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

