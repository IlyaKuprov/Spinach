% SABRE experiment simulation for Eibe Duecker and Christian Griesinger.
% Set to reproduce Figure 3b from http://dx.doi.org/10.1021/ja903601p
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function sabre_pyridine()

% Spin system
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H'};

% Chemical shifts
inter.zeeman.scalar={8.54 7.44 7.86 7.44 8.54 -23.5 -23.5};

% Couplings inside pyridine
inter.coupling.scalar=cell(7);
inter.coupling.scalar{1,2}= 4.88; 
inter.coupling.scalar{4,5}= 4.88;
inter.coupling.scalar{1,4}= 1.00; 
inter.coupling.scalar{2,5}= 1.00;
inter.coupling.scalar{1,3}= 1.84; 
inter.coupling.scalar{3,5}= 1.84;
inter.coupling.scalar{1,5}=-0.13;
inter.coupling.scalar{2,3}= 7.67; 
inter.coupling.scalar{3,4}= 7.67;
inter.coupling.scalar{2,4}= 1.37;

% Couplings of the hydride group 
inter.coupling.scalar{1,6}=1.12;
inter.coupling.scalar{1,7}=1.02;
inter.coupling.scalar{6,7}=7.00;

% Magnetic fields
polarization_field=25e-3;
nmr_field=7.05;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.enable={'greedy'};

% Do the housekeeping
sys.magnet=polarization_field;
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get the Hamiltonian superoperator
H=hamiltonian(assume(spin_system,'nmr'));

% Start in a singlet state
rho=singlet(spin_system,6,7);

% Evolve the system for 2.5 seconds
rho=evolution(spin_system,H,[],rho,2.5,1,'final');

% Disconnect the parahydrogen
[H,rho]=decouple(spin_system,H,rho,[6 7]);

% Evolve the system for a further 2.5 seconds
rho=evolution(spin_system,H,[],rho,2.5,1,'final');

% Split the Hamiltonian superoperator
Hz=hamiltonian(assume(spin_system,'nmr','zeeman'));
Hc=hamiltonian(assume(spin_system,'nmr','couplings'));
Hz=decouple(spin_system,Hz,[],[6 7])/polarization_field;
Hc=decouple(spin_system,Hc,[],[6 7]);

% Lift the field exponentially in 1024 steps over 5 seconds
 for field=exp(linspace(log(polarization_field),log(nmr_field),1024))
     rho=step(spin_system,Hc+field*Hz,rho,5/1024);
 end

% Set the field to NMR field
H=Hc+nmr_field*Hz;
 
% Evolve the system for a further second in high field
rho=evolution(spin_system,H,[],rho,1.0,1,'final');

% Set NMR experiment parameters
parameters.offset=2400;
parameters.spins={'1H'};
parameters.sweep=600;
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.axis_units='kHz';
parameters.invert_axis=1;

% Set the detection state
coil=state(spin_system,'L+',parameters.spins{1});

% Get the pulse operator
Ly=operator(spin_system,'Ly',parameters.spins{1});

% Apply pi/2 pulse on Y axis
rho=step(spin_system,Ly,rho,pi/2);

% Detect the magnetization
fid=evolution(spin_system,H,coil,rho,1/parameters.sweep,...
              parameters.npoints-1,'observable');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,-real(spectrum),parameters);

end

