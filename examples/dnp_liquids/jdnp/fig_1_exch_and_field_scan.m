% Matching condition plot for JDNP - proton polarisation at
% a particular time as a function of the external field and
% the inter-electron exchange coupling. Further details in:
%
%           https://doi.org/10.1039/d1cp04186j
%
% Calculation time: hours, line-by-line plotting 
%
% maria-grazia.concilio@weizmann.ac.il
% i.kuprov@soton.ac.uk

function fig_1_exch_and_field_scan()

% Load the spin system
[sys,inter,bas,parameters]=system_specification();

% Experiment parameters
parameters.mw_pwr=(2*pi*1e6)/2;       % rad/s
parameters.mw_dur=20e-3;              % seconds

% Field and coupling grids
field_grid=linspace(0.25,3.0,256);    % Tesla
exch_grid=linspace(-100e9,100e9,256); % Hz

% Preallocate output array
dnp=zeros([numel(field_grid) numel(exch_grid)]);

% Create and scale the figure
current_fig=figure(); scale_figure([1.0 0.65]);

% Loop over the fields
for n=1:numel(field_grid)

    % Set the magnet field
    sys.magnet=field_grid(n);  
                                               
    % Trityl and free electron frequencies 
    f_free=2*pi*g2freq(parameters.g_ref,sys.magnet);
    f_trityl=2*pi*g2freq(parameters.g_trityl,sys.magnet);
    
    % Microwave offset
    parameters.mw_off=f_trityl-f_free;      
            
    % Parallel loop over exchange couplings
    parfor k=1:numel(exch_grid) %#ok<*PFBNS>

        % Localisation
        locinter=inter;

        % Set the exchange coupling
        locinter.coupling.scalar=cell(3,3); 
        locinter.coupling.scalar{2,3}=exch_grid(k);
  
        % Spinach housekeeping
        spin_system=create(sys,locinter);  
        spin_system=basis(spin_system,bas);
    
        % Relevant operators and states
        E_mw=operator(spin_system,'Lx','E')/2;
        Ez=operator(spin_system,'Lz','E');
      
        % Thermal equilibium state
        H0=hamiltonian(assume(spin_system,'labframe'),'left');
        rho_eq=equilibrium(spin_system,H0); 

        % Detection state
        Nz=state(spin_system,'Lz','1H');

        % Hamiltonian and relaxation superoperator
        H=hamiltonian(assume(spin_system,'esr'));
        R=relaxation(spin_system);
    
        % Add microwave operator
        H=H+parameters.mw_pwr*E_mw;  

        % Add microwave offset
        H=H+parameters.mw_off*Ez;
        
        % Run the time evolution
        rho=evolution(spin_system,H+1i*R,[],rho_eq,...
                      parameters.mw_dur,1,'final');
                       
        % Compute the enhancement vs equilibrium
        dnp(n,k)=real(Nz'*rho)/real(Nz'*rho_eq);
        
    end
    
    % Update the plot
    set(groot,'CurrentFigure',current_fig);
    imagesc(exch_grid/1e9,field_grid,dnp);
    kxlabel('$\omega_{ex}$ / GHz'); kylabel('$B_0$, Tesla');
    kcolourbar('$^{1}$H DNP @ 20 ms'); drawnow;
        
end

end

