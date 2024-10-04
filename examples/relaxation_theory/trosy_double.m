% Hari Arthanari's Double TROSY effect.
% 
% Calculation time: seconds.
%
% i.kuprov@soton.ac.uk

function trosy_double()

% Magnet field
sys.magnet=14.1;

% Read 3-fluorotyrosine DFT calculation
[~,inter_dft]=g2spinach(gparse('../standard_systems/4_fluoro_phe.out'),...
              {{'H','1H'},{'C','13C'},{'F','19F'}},[32.07 186.38 192.97]);
                                 
% Extract coordinates and CSAs
sys.isotopes={'1H','1H','19F','13C'};
inter.zeeman.matrix=cell(1,4);
inter.zeeman.matrix{1}=inter_dft.zeeman.matrix{10};
inter.zeeman.matrix{2}=inter_dft.zeeman.matrix{20};
inter.zeeman.matrix{3}=inter_dft.zeeman.matrix{19};
inter.zeeman.matrix{4}=inter_dft.zeeman.matrix{8};
inter.coordinates=cell(4,1);
inter.coordinates{1}=inter_dft.coordinates{10};
inter.coordinates{2}=inter_dft.coordinates{20};
inter.coordinates{3}=inter_dft.coordinates{19};
inter.coordinates{4}=inter_dft.coordinates{8};

% Extract J-couplings
inter.coupling.scalar=inter_dft.coupling.scalar([10 20 19 8],[10 20 19 8]);

% Relaxation theory
inter.relaxation={'redfield'};
inter.rlx_keep='secular';
inter.equilibrium='zero';
inter.tau_c={20e-9};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters - 13C
parameters.spins={'13C'};
parameters.rho0=state(spin_system,'L+','13C');
parameters.coil=state(spin_system,'L+','13C');
parameters.decouple={};
parameters.offset=26800;
parameters.sweep=500;
parameters.npoints=2048;
parameters.zerofill=16384;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'gauss',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

                         