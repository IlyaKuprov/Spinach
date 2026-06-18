% INADEQUATE spectrum of a three-spin system with J-coupling between
% two spins only. The sequence selects double-quantum coherence
% from coupled 13C pairs and converts it back for detection.
%
% Calculation time: seconds
% 
% bud.macaulay@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function inad_three_spin()

% Magnet field
sys.magnet=9.4;

% Spin system and interactions
sys.isotopes={'13C','13C','13C'};
inter.zeeman.scalar={10 15 20};
inter.coupling.scalar{1,2}=55;
inter.coupling.scalar{3,3}=0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'13C'};
parameters.J=55;
parameters.decouple={};
parameters.offset=1800;
parameters.sweep=5000;
parameters.npoints=4096;
parameters.zerofill=16384;
parameters.axis_units='ppm';

% Simulation
fid=liquid(spin_system,@inadequate,parameters,'nmr');

% Processing
fid=apodisation(spin_system,fid,{{'exp',5}});
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

