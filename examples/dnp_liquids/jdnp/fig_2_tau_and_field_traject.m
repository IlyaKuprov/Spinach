% Time evolution plot for JDNP: proton polarisation as a function 
% of time for specific external fields and rotational correlation 
% times. The inter-electron exchange coupling is set to the match-
% ing condition at each field. Further details in:
%
%               https://doi.org/10.1039/d1cp04186j
%
% Calculation time: seconds, line-by-line plotting 
%
% maria-grazia.concilio@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il

function fig_2_tau_and_field_traject()

% Load the spin system
[sys,inter,bas,parameters]=system_specification();

% Experiment parameters
parameters.mw_pwr=2*pi*250e3;
parameters.t_step=1e-3;
parameters.nsteps=200;    

% Magnetic field grid, Tesla
field_grid=[0.5 3.4 7.0 11.7 14.1 23.5];

% Correlation time grid, seconds
tau_c=[300e-12 400e-12 500e-12 600e-12];

% Get a figure going
figure(); scale_figure([1.75 1.0]);

% Loop over the field grid
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
    
    % Set figure properties
    subplot(2,3,n); hold on;
    axis tight; ylim([-250 50]); kgrid; box on;
    ktitle(['$B_0$ = ' num2str(field_grid(n)) ' Tesla']);
    kylabel('$^{1}$H DNP'); kxlabel('time / ms');
    
    % Loop over correlation times
    for k=1:numel(tau_c)
      
        % Set correlation time      
        inter.tau_c={tau_c(k)};
  
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
        plot(t_axis,answer); drawnow;
             
    end
    
    % Set the legend
    klegend({'$\tau_c = 300$ ps','$\tau_c = 400$ ps',...
             '$\tau_c = 500$ ps','$\tau_c = 600$ ps'},...
             'location','southeast');
    
end

end

