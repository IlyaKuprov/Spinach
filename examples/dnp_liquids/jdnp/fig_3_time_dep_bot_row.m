% Time evolution plot for JDNP: proton polarisation as a function 
% of time for specific external fields. The inter-electron exchan-
% ge coupling is set to the matching condition at each field. See
% also the "top row" simulation where one of the electrons is re-
% moved to demostrate that JDNP vanishes. Further details in:
%
%               https://doi.org/10.1039/d1cp04186j
%
% Calculation time: seconds, line-by-line plotting 
%
% maria-grazia.concilio@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il

function fig_3_time_dep_bot_row()

% Load the spin system
[sys,inter,bas,parameters]=system_specification();

% Experiment parameters
parameters.mw_pwr=2*pi*250e3;
parameters.t_step=1e-3;
parameters.nsteps=300;    

% Magnetic field grid, Tesla
field_grid=[0.034 0.34 3.4];

% Get a figure going
figure(); scale_figure([1.75 0.5]);

% Loop over the fields
for n=1:numel(field_grid)
   
    % Set magnet field
    sys.magnet=field_grid(n);
    
    % Set microwave offset frequency
    f_free=g2freq(parameters.g_ref,sys.magnet);
    f_trityl=g2freq(parameters.g_trityl,sys.magnet);
    parameters.mw_off=2*pi*(f_trityl-f_free);      
    
    % Set the exchange coupling
    electron_zeeman_tensor=g2freq(inter.zeeman.matrix{2},sys.magnet);
    electron_zeeman_iso=mean(diag(electron_zeeman_tensor));
    proton_zeeman_iso=sys.magnet*spin('1H')/(2*pi);
    inter.coupling.scalar{2,3}=electron_zeeman_iso+proton_zeeman_iso;
        
    % Spinach housekeeping
    spin_system=create(sys,inter);  
    spin_system=basis(spin_system,bas);
    
    % Electron operators
    Ex=operator(spin_system,'Lx','E');
    Ez=operator(spin_system,'Lz','E');
        
    % Thermal equilibirium state
    H0=hamiltonian(assume(spin_system,'labframe'),'left');
    rho_eq=equilibrium(spin_system,H0);
        
    % Hamiltonian and relaxation superoperator
    H=hamiltonian(assume(spin_system,'esr'));
    R=relaxation(spin_system);
    
    % Add microwave terms to the Hamiltonian
    H=H+parameters.mw_pwr*Ex+parameters.mw_off*Ez;
             
    % Detection state: Lz on the proton
    coil=state(spin_system,'Lz','1H');
             
    % Run the time evolution and normalise to thermal equilibrium
    answer=evolution(spin_system,H+1i*R,coil,rho_eq,parameters.t_step,...
                     parameters.nsteps,'observable');
    answer=real(answer/(coil'*rho_eq));
   
    % Time axis generation
    t_axis=linspace(0,parameters.t_step*parameters.nsteps,...
                    parameters.nsteps+1);  
        
    % Plotting
    subplot(1,3,n); plot(t_axis,answer); kgrid; box on;
    ktitle(['$B_0$ = ' num2str(field_grid(n)) ' Tesla']);
    kylabel('$^{1}$H DNP'); kxlabel('time / seconds'); drawnow;
  
end

end

