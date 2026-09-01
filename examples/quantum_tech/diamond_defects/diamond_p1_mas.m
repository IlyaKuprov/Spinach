% Powder magic angle spinning EPR spectrum of the P1 substitutional
% nitrogen defect in diamond at 7 Tesla. The initial condition is a
% transverse electron magnetisation, and therefore no pulse is need-
% ed: the single angle spinning context calls acquire.m directly.
%
% The anisotropy is carried by the 14N hyperfine coupling, of which
% the dipolar part is 10.9 MHz; the two outer hyperfine lines there-
% fore acquire sideband manifolds, whereas the central line has no
% hyperfine anisotropy and stays sharp. The rotation frequency used
% here is above what a mechanical rotor can deliver, but within the
% range reported for optically levitated nanoparticles:
%
%           R. Reimann et al., Phys. Rev. Lett. 121, 033602 (2018)
%           J. Ahn et al., Phys. Rev. Lett. 121, 033603 (2018)
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function diamond_p1_mas()

% P1 centre parameters
p1_params.orientation='111';
p1_params.nitrogen='14N';

% Build the spin system
[sys,inter]=diamond_p1(p1_params);

% Magnet field
sys.magnet=7;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Rotor parameters
parameters.rate=5e6;
parameters.axis=[1 1 1];
parameters.max_rank=17;

% Sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'L+','E');
parameters.coil=state(spin_system,'L+','E');
parameters.offset=1.234e7;
parameters.sweep=3e8;
parameters.npoints=4096;
parameters.zerofill=16384;
parameters.grid='rep_2ang_1600pts_sph';
parameters.axis_units='MHz-labframe';
parameters.invert_axis=1;
parameters.verbose=0;

% Simulation
fid=singlerot(spin_system,@acquire,parameters,'esr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

