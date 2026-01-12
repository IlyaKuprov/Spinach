% Static powder 79Br NMR spectrum of potassium bromide. At least
% 3 quadrupolar tensors are necessary to reproduce the experimen-
% tal shape, likely due to a distribution of electrostatic envi-
% ronments in the powder.
%
% Calculation time: seconds.
%
% sanjay.vinod-kumar@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il

function static_powder_nqi_b()

% Magnet field
sys.magnet=9.3659;

% Spin system
sys.isotopes={'79Br','79Br','79Br'};

% Chemical shift, ppm
inter.zeeman.scalar={60.0933 60.0933 60.0933};

% Quadrupolar coupling
inter.coupling.matrix{1,1}=1e3*diag([13.7569 1.6424 -(13.7569 + 1.6424)]);
inter.coupling.matrix{2,2}=1e3*diag([4.0779  4.5179 -( 4.0779 + 4.5179)]);
inter.coupling.matrix{3,3}=1e3*diag([1.5885  0.9449 -( 1.5885 + 0.9449)]);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=1; bas.projections=+1;

% Algorithmic options
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.grid='icos_2ang_163842pts';
parameters.spins={'79Br'};
parameters.decouple={};
parameters.axis_units='Hz';
parameters.invert_axis=1;
parameters.offset=6034.96;  % receiver offset, Hz
parameters.sweep=1e5;       % sweep width, Hz
parameters.npoints=1024;     % points to acquire
parameters.zerofill=4096;   % zerofill to
parameters.rho0=40*state(spin_system,{'L+'},{1})+...
                32*state(spin_system,{'L+'},{2})+...
                28*state(spin_system,{'L+'},{3});
parameters.coil=state(spin_system,'L+','79Br');

% Simulation
fid=powder(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spec=real(fftshift(fft(fid,parameters.zerofill)));

% Plotting
kfigure(); plot_1d(spin_system,spec,parameters);
ylim([-10 1000]);

end

