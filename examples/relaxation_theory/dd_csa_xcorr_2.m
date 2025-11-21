% DD-CSA cross-correlation - a reproduction of Fig 5a from the paper
% by Grace and Kumar (http://dx.doi.org/10.1006/jmra.1995.1151).
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function dd_csa_xcorr_2()

% Read the spin system parameters (vacuum DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/fdnb.log'),...
                   {{'H','1H'},{'F','19F'}},[32.0 270.0],[]);

% Set up the calculation
sys.magnet=9.4;                  % Magnet induction
sys.tols.prox_cutoff=5;          % Increase proximity cutoff to 5 Angstrom
inter.relaxation={'redfield'};   % Redfield relaxation theory
inter.rlx_keep='secular';        % Rotating-frame version of Redfield theory
inter.equilibrium='dibari';      % Relax to thermal equilibrium
inter.temperature=298;           % Work at room tempearture
inter.tau_c={9.6e-12};           % Correlation time

bas.formalism='sphten-liouv';    % Liouville space formalism
bas.approximation='none';        % Complete basis set

% Proximity cut-off
sys.tols.prox_cutoff=4.0;

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set simulation parameters
parameters.offset=-521;       % Axis offset
parameters.sweep=50;          % Sweep width
parameters.npoints=128;       % Number of simulation points
parameters.zerofill=512;      % Zerofill fid to twice the length 
parameters.spins={'19F'};     % Run the experiment on 19F

% Set the assumptions to high-field NMR
spin_system=assume(spin_system,'nmr');

% Get the Hamiltonian superoperator
L=hamiltonian(spin_system);

% Add Redfield superoperator,
L=L+1i*relaxation(spin_system);

% Apply the offset
L=frqoffset(spin_system,L,parameters);

% Get pulse operators on fluorine
Lx=operator(spin_system,'Lx','19F');
Ly=operator(spin_system,'Ly','19F');

% Set detection state to L+
coil=state(spin_system,'L+','19F');

% Calculate the time step of the simulation
timestep=1/parameters.sweep;

% Compute the axis
axis_hz=sweep2ticks(parameters.offset,parameters.sweep,parameters.zerofill);

% Set the initial state to thermal equilibrium
rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'));

% Get the figure going
kfigure(); 

% Loop over mixing times
for t_mix=[0.1 1.4 1.6 1.8 2.0 2.2 2.4 10]
    
    % Apply the inversion pulse
    rho=step(spin_system,Lx,rho_eq,pi);
    
    % Run the mixing time
    rho=evolution(spin_system,L,coil,rho,t_mix,1,'final');
    
    % Apply the detection pulse
    rho=step(spin_system,Ly,rho,pi/2);
    
    % Run the detection period
    fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable');
    
    % Perform the Fourier transform
    spectrum=fftshift(fft(apodisation(spin_system,fid,{{'exp',6}}),parameters.zerofill));
    
    % Plot the spectrum
    hold on; plot(axis_hz,real(spectrum)); drawnow;
    
end

% Invert the X axis
set(gca,'XDir','reverse'); 
xlim tight; ylim padded; box on; 
kgrid; kxlabel('19F linear frequency, Hz');
klegend('$\tau_m = 0.1$ s','$\tau_m = 1.4$ s','$\tau_m = 1.6$ s',...
        '$\tau_m = 1.8$ s','$\tau_m = 2.0$ s','$\tau_m = 2.2$ s',...
        '$\tau_m = 2.4$ s','$\tau_m = 10$ s');

end

