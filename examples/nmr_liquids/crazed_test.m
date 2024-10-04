% Long range intermolecular coherences predicted by Warren and 
% co-workers
%
%         http://dx.doi.org/10.1126/science.8266096
%
% Calculation time: seconds
%
% hannah.hogben@chem.ox.ac.uk
% i.kuprov@soton.ac.uk

function crazed_test()

% Specify system parameters
sys.magnet=6.0;
sys.isotopes={'1H','1H','1H','1H'};
inter.zeeman.scalar={2.0,2.0,8.0,8.0};
inter.coordinates={[0 0 0];
                   [1/2 1/2 -1/sqrt(2)]*1e2;
                   [1 0 0]*1e2;
                   [1/2 -1/2 -1/sqrt(2)]*1e2};
inter.temperature=100;
                              
% Use the complete basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach code
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.angle=pi/2;
parameters.offset=1300;
parameters.sweep=5000;
parameters.npoints=[512 512];
parameters.zerofill=[2048 2048];
parameters.spins={'1H'};
parameters.orientation=[0 0 0];

% Thermal equilibrium state
[H,Q]=hamiltonian(assume(spin_system,'labframe'),'left');
parameters.rho0=equilibrium(spin_system,H,Q,[0 0 0]);

% CRAZED simulation
fid=crystal(spin_system,@crazed,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'cos'},{'cos'}});

% Fourier transform
spectrum=fftshift(fft2(fid,parameters.zerofill(2),...
                           parameters.zerofill(1)));

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,abs(spectrum),parameters,...
        20,[0.02 0.2 0.02 0.2],2,256,6,'positive');

end

