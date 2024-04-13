% Optimisation of slice durations in a composite inversion
% pulse with specified amplitudes, phases, and a constrain-
% ed overall duration. The initial guess is
%
%      270(-x)360(x)90(y)270(-y)360(y)90(x)  [Fig. 3]
%
% from https://doi.org/10.1016/0022-2364(83)90133-6 -- the
% optimisation demonstrates that a slightly better pulse of
% the same power and duration exists.
%
% Calculation time: minutes.
%
% i.kuprov@soton.ac.uk

function features_dt_opt()

% Set the magnetic field
sys.magnet=14.1;

% Put 100 non-interacting spins at equal intervals over the area
% that needs to be affected by the pulse (25 kHz either side)
n_spins=100; sys.isotopes=cell(n_spins,1);
for n=1:n_spins
    sys.isotopes{n}='13C';
end
inter.zeeman.scalar=num2cell(linspace(-166,166,n_spins));

% Select a basis set - IK-2 keeps complete basis on each 
% spin in this case, but ignores multi-spin orders
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.space_level=1;
bas.connectivity='scalar_couplings';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up initial and target states
rho_init=state(spin_system,'Lz','13C'); 
rho_init=rho_init/norm(full(rho_init),2);
rho_targ=rho_init; % Inversion when minimised

% Get the control operators
controls{1}=operator(spin_system,'Lx','13C');
controls{2}=operator(spin_system,'Ly','13C');

% Get the drift Hamiltonian
drift=hamiltonian(assume(spin_system,'nmr'));

% Specify the pulse sequence (rad/s)
waveform=2*pi*2.5e4*[-1  1  0  0  0  1;
                      0  0  1 -1  1  0];
dt_old=10*[3; 4; 1; 3; 4; 1]; time_unit=1e-6;

% Optimisation options
options=optimoptions('fmincon','SpecifyObjectiveGradient',true,...
                     'Display','iter','HessianApproximation','lbfgs');

% Use Matlab's optimiser with a total duration constraint
target_fun=@(x)tgrape(spin_system,drift,controls,waveform,...
                      x,time_unit,rho_init,rho_targ);
lb=zeros(size(dt_old)); Aeq=ones(size(dt_old))'; beq=sum(dt_old);
dt_new=fmincon(target_fun,dt_old,[],[],Aeq,beq,lb,[],[],options);

% Apply the old and the new pulse
rho_old=rho_init; rho_new=rho_init;
Lx=controls{1}; Ly=controls{2};
for n=1:numel(dt_new)
    rho_old=step(spin_system,drift+waveform(1,n)*Lx+waveform(2,n)*Ly,...
                 rho_old,dt_old(n)*time_unit);
    rho_new=step(spin_system,drift+waveform(1,n)*Lx+waveform(2,n)*Ly,...
                 rho_new,dt_new(n)*time_unit);
end

% Set acquisition parameters
parameters.spins={'13C'};
parameters.coil=state(spin_system,'L+','13C');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=55000;
parameters.npoints=2048;
parameters.zerofill=16384;
parameters.axis_units='Hz';
parameters.invert_axis=1;

% Run the simulation for the initial guess
parameters.rho0=step(spin_system,Ly,rho_old,pi/2);
fid=liquid(spin_system,@acquire,parameters,'nmr');
fid=apodization(fid,'gaussian-1d',10);
spectrum_old=fftshift(fft(fid,parameters.zerofill));

% Run the simulation for the optimised pulse
parameters.rho0=step(spin_system,Ly,rho_new,pi/2);
fid=liquid(spin_system,@acquire,parameters,'nmr');
fid=apodization(fid,'gaussian-1d',10);
spectrum_new=fftshift(fft(fid,parameters.zerofill));

% Plot the spectrum before and after optimisation
figure(); scale_figure([2.0 0.75]);
subplot(1,2,1); plot_1d(spin_system,real(spectrum_old),parameters);
ktitle('composite inversion, initial guess');
subplot(1,2,2); plot_1d(spin_system,real(spectrum_new),parameters);
ktitle('composite inversion, optimised');

% Write the old and the new pulses to the console
disp('Old and new durations (s):'); disp([dt_old dt_new]);

end

