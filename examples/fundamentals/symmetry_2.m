% Pulse-acquire NMR spectrum of a highly symmetric spin 
% system provided by Andres Castillo. Uses the fully sym-
% metric irreducible representation of S3(x)S3(x)S3 per-
% mutation symmetry group.
%
% i.kuprov@soton.ac.uk

function symmetry_2()

% Spin system specification
sys.magnet=9.4;
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'};

inter.zeeman.scalar={0.89 0.89 0.89 0.895 0.895 0.895 1.16 1.16 1.16 1.2 1.39 1.71 3.85};
inter.coupling.scalar=num2cell(...
      [  0         0         0         0         0         0         0         0         0         0         0    6.7200         0
         0         0         0         0         0         0         0         0         0         0         0    6.7200         0
         0         0         0         0         0         0         0         0         0         0         0    6.7200         0
         0         0         0         0         0         0         0         0         0         0         0    6.6400         0
         0         0         0         0         0         0         0         0         0         0         0    6.6400         0
         0         0         0         0         0         0         0         0         0         0         0    6.6400         0
         0         0         0         0         0         0         0         0         0         0         0         0    6.0800
         0         0         0         0         0         0         0         0         0         0         0         0    6.0800
         0         0         0         0         0         0         0         0         0         0         0         0    6.0800
         0         0         0         0         0         0         0         0         0         0   14.0850    5.1000    8.3000
         0         0         0         0         0         0         0         0         0   14.0850         0    8.3050    6.4000
    6.7200    6.7200    6.7200    6.6400    6.6400    6.6400         0         0         0    5.1000    8.3050         0         0
         0         0         0         0         0         0    6.0800    6.0800    6.0800    8.3000    6.4000         0         0]/2);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.sym_group={'S3','S3','S3'};
bas.sym_spins={[1 2 3],[4 5 6],[7 8 9]};
bas.connectivity='scalar_couplings';
bas.space_level=1;
bas.projections=+1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
     
% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=800;
parameters.sweep=2000;
parameters.npoints=2048;
parameters.zerofill=8196;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);     

end

