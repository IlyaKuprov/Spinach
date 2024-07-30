% 1H NMR spectrum of a large and highly symmetric spin system
% with two tert-butyl groups supplied by Eberhard Matern. Done
% by brute force in Hilbert space.
%
% WARNING: needs 32 CPU cores, 128 GB of RAM and
%          a Titan V or later.
%
% Run time on the above: hours
%
% i.kuprov@soton.ac.uk
% eberhard.matern@t-online.de

function high_symmetry_1()

% Isotopes
sys.isotopes={'31P','31P','1H','1H',...
              '1H','1H','1H','1H','1H','1H','1H','1H','1H',...
              '1H','1H','1H','1H','1H','1H','1H','1H','1H'};

% Magnetic induction
sys.magnet=9.39798;

% Chemical shifts
inter.zeeman.scalar={-43.844, -43.844, 4.090, 4.090,...
                       1.354,   1.354, 1.354, 1.354, 1.354, 1.354 1.354, 1.354, 1.354,...
                       1.354,   1.354, 1.354, 1.354, 1.354, 1.354 1.354, 1.354, 1.354};

% Scalar couplings
inter.coupling.scalar=cell(22);
inter.coupling.scalar{1,2}=301.99; 
inter.coupling.scalar{1,3}=-321.62; 
inter.coupling.scalar{2,4}=-321.62; 
inter.coupling.scalar{1,4}=-19.15; 
inter.coupling.scalar{2,3}=-19.15; 

inter.coupling.scalar{1,5}=15.63; 
inter.coupling.scalar{1,6}=15.63; 
inter.coupling.scalar{1,7}=15.63; 
inter.coupling.scalar{1,8}=15.63; 
inter.coupling.scalar{1,9}=15.63; 
inter.coupling.scalar{1,10}=15.63; 
inter.coupling.scalar{1,11}=15.63; 
inter.coupling.scalar{1,12}=15.63; 
inter.coupling.scalar{1,13}=15.63; 

inter.coupling.scalar{2,14}=15.63; 
inter.coupling.scalar{2,15}=15.63; 
inter.coupling.scalar{2,16}=15.63; 
inter.coupling.scalar{2,17}=15.63; 
inter.coupling.scalar{2,18}=15.63; 
inter.coupling.scalar{2,19}=15.63; 
inter.coupling.scalar{2,20}=15.63; 
inter.coupling.scalar{2,21}=15.63; 
inter.coupling.scalar{2,22}=15.63; 

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Algorithmic options
sys.enable={'gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=1150;
parameters.sweep=2400;
parameters.npoints=4096;
parameters.zerofill=32768;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodization
fid=apodization(fid,'exp-1d',5);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

