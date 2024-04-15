% An illustration of basis set coefficient optimisation for a pul-
% se optimised as a linear combination of user-specified vectors.
% The optimisation is done using LBFGS GRAPE algorithm.
%
%          http://dx.doi.org/10.1016/j.jmr.2011.07.023
%
% In a set of 100 equspaced signals, the central 60 spins are set
% up for maximum excitation; there are no constraints on the dyna-
% mics of the 20 spins on either side of the interval.
%
% Calculation time: minutes.
%
% i.kuprov@soton.ac.uk
% david.goodwin@inano.au.dk

function features_wave_basis()

% Set the magnetic field
sys.magnet=14.1;

% 100 non-interacting spins at equal intervals 
% within the plus/minus 160 ppm range
n_spins=100; sys.isotopes=cell(n_spins,1);
for n=1:n_spins, sys.isotopes{n}='13C'; end
inter.zeeman.scalar=num2cell(linspace(-160,160,n_spins));

% Select a basis set - IK-2 keeps complete basis on each 
% spin in this case, but ignores multi-spin orders
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.space_level=1;
bas.connectivity='scalar_couplings';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up the initial state vector
rho_init=state(spin_system,'Lz',21:80);    % 60 spins in the centre
rho_init=rho_init/norm(full(rho_init),2);

% Set up the target state vector
rho_targ=state(spin_system,'Lx',21:80);    % 60 spins in the centre
rho_targ=rho_targ/norm(full(rho_targ),2);

% Get the control operators
Cx=operator(spin_system,'Lx','13C');
Cy=operator(spin_system,'Ly','13C');

% Get the drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Define control parameters
control.drifts={{H}};                            % Drift
control.operators={Cx,Cy};                       % Controls
control.rho_init={rho_init};                     % Starting state
control.rho_targ={rho_targ};                     % Target state
control.pulse_dt=4e-6*ones(1,125);               % Slice durations
control.basis=wave_basis('legendre',20,125)';    % Basis set
control.pwr_levels=2*pi*linspace(15e3,20e3,11);  % Power levels
control.method='lbfgs';                          % Optimisation method
control.max_iter=200;                            % Termination condition
control.penalties={'SNS'};                       % Penalty type
control.p_weights=10;                            % Penalty weight
control.parallel='ensemble';                     % Parallelisation 

% Plotting options
control.plotting={'xy_controls','amp_controls',...
                  'robustness','spectrogram'};

% Initial guess
guess=randn(2,20)/20;

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Run LBFGS GRAPE pulse optimisation
basis_coeffs=fminnewton(spin_system,@grape_xy,guess);

% Reassemble time-domain control sequence
pulse=mean(control.pwr_levels)*basis_coeffs*control.basis;
pulse=mat2cell(pulse,[1 1]);

% Simulate the optimised pulse
rho_init=state(spin_system,'Lz','13C');
rho=shaped_pulse_xy(spin_system,H,{Cx,Cy},pulse,...
                    control.pulse_dt,rho_init,'expv-pwc');

% Set acquisition parameters
parameters.spins={'13C'};
parameters.rho0=rho;
parameters.coil=state(spin_system,'L+','13C');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=55000;
parameters.npoints=2048;
parameters.zerofill=16384;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulate the free induction decay
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodization
fid=apodization(fid,'gaussian-1d',10);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);
hold on; xline(-96,'r-'); xline(96,'r-');
    
end

