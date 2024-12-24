% A demonstration that the JDNP effect vanishes when the second 
% electron is removed from the system. Proton polarisation as a
% function of time for specific external fields is plotted. See
% also the "bot row" simulation where both electrons are active
% and the JDNP enhancement is present. Further details in:
%
%               https://doi.org/10.1039/d1cp04186j
%
% Calculation time: seconds, line-by-line plotting 
%
% maria-grazia.concilio@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il

function fig_3_time_dep_top_row()

% Load the spin system
[sys,inter,bas,parameters]=system_specification();

% Kill the second electron 
sys.isotopes=sys.isotopes([1 2]);
inter.zeeman.matrix=inter.zeeman.matrix([1 2]);
inter.coordinates=inter.coordinates([1 2]);
inter.coupling.scalar=cell(2,2);
inter=rmfield(inter,{'srfk_tau_c','srfk_mdepth'});
inter.relaxation={'redfield'};

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

