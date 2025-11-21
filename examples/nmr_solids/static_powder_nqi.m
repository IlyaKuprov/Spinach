% Static quadrupolar 14N powder pattern of L-valyl-L-alanine using 
% very large numerical orientation grid, set to reproduce Figure 5
% from the paper by O'Dell and Ratcliffe:
%
%           http://dx.doi.org/10.1016/j.cplett.2011.08.030
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function static_powder_nqi()

% System specification
sys.magnet=21.1;
sys.isotopes={'14N','14N'};
inter.coupling.matrix{1,1}=eeqq2nqi(1.24e6,0.22,1,[0 0 0]);
inter.coupling.matrix{2,2}=eeqq2nqi(3.06e6,0.40,1,[0 0 0]);

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'14N'};
parameters.decouple={};
parameters.offset=0;
parameters.sweep=6e6;
parameters.npoints=512;
parameters.zerofill=2048;
parameters.axis_units='MHz';
parameters.invert_axis=1;
parameters.grid='icos_2ang_163842pts';
parameters.rho0=state(spin_system,'L+','14N');
parameters.coil=state(spin_system,'L+','14N');
parameters.verbose=0;

% Simulation
fid=powder(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

