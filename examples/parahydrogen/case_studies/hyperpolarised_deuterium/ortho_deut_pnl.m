% Simulated PNL (partially negative line) spectrum of 
% ortho-deuterium in the presence of a parahydrogena-
% tion catalyst. Paper link coming in due course.
%
% thhu@mpinat.mpg.de
% anakin.aden@mpinat.mpg.de
% denismoll@hotmail.de
% julius.matz@mpinat.mpg.de
% stefan.gloeggler@mpinat.mpg.de
% ilya.kuprov@weizmann.ac.il

function ortho_deut_pnl()

% Magnet field
sys.magnet=7.05;

% Spin system
sys.isotopes={'2H','2H','2H','2H'};
inter.zeeman.scalar={4.55 4.55 -13.5 -16.5};
inter.coupling.scalar=cell(4,4);
inter.coupling.scalar{1,2}=12.0; % free deuterium
inter.coupling.scalar{3,4}=0.25; % bound to catalyst

% Kinetics
inter.chem.parts={[1 2],[3 4]};
inter.chem.rates=[-1  7000; 
                   1 -7000];
inter.chem.concs=[1 0];

% Simulation formalsim
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.sweep=200;
parameters.offset=200;
parameters.npoints=1024;
parameters.spins={'2H'};
parameters.axis_units='ppm';
parameters.zerofill=4096;

% Initial condition
options.dephasing=1;
[S,~,Q]=deut_pair(spin_system,1,2,options);
rho0=S+Q{1}+Q{2}+Q{3}+Q{4}+Q{5};

% Detection state
coil=state(spin_system,'L+','2H'); 

% Build the Liouvillian
spin_system=assume(spin_system,'nmr');
L=hamiltonian(spin_system)+1i*kinetics(spin_system);
L=frqoffset(spin_system,L,parameters);

% Build the pulse operator
Dy=operator(spin_system,'Ly','2H');

% Apply the 45-degree pulse
rho=step(spin_system,Dy,rho0,pi/4);

% Run the evolution and watch the coil state
fid=evolution(spin_system,L,coil,rho,1/parameters.sweep,...
              parameters.npoints-1,'observable');

% Apodisation and Fourier transform
fid=apodisation(spin_system,fid,{{'exp',6}});
spec=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spec),parameters);
xlim([4.4 4.7]); ylim padded; kylabel('intensity');
kxlabel('$^2$H chemical shift / ppm'); kgrid;

end

