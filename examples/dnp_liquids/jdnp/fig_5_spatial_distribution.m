% An illustration of the fact that JDNP effect does not vanish on
% position and orientation averaging in liquid phase. The simula-
% tion shows proton polarisation at 20 ms for different proton lo-
% cations around a radical pair with inter-electron exchange cou-
% pling chosen to acheive the JDNP effect. Details in:
%
%               https://doi.org/10.1039/d1cp04186j
%
% Calculation time: seconds
%
% mariagrazia.concilio@sjtu.edu.cn
% ilya.kuprov@weizmann.ac.il

function fig_5_spatial_distribution()

% Load the spin system
[sys,inter,bas,parameters]=system_specification();

% Set magnet field
sys.magnet=14.08;

% Experiment parameters
parameters.mw_pwr=2*pi*250e3;
parameters.t_evol=20e-3;

% Set microwave offset frequency
f_free=g2freq(parameters.g_ref,sys.magnet);
f_trityl=g2freq(parameters.g_trityl,sys.magnet);
parameters.mw_off=2*pi*(f_trityl-f_free);   

% Set the exchange coupling
electron_zeeman_tensor=g2freq(inter.zeeman.matrix{2},sys.magnet);
electron_zeeman_iso=mean(diag(electron_zeeman_tensor));
proton_zeeman_iso=sys.magnet*spin('1H')/(2*pi);
inter.coupling.scalar{2,3}=electron_zeeman_iso+proton_zeeman_iso;

% Specify coordinate arrays
X=linspace(-30,30,30); Y=linspace(-30,30,30); Z=linspace(-30,30,30);

% Get a figure going
figure(); scale_figure([1.75 1.0]);

% Preallocate answer arrays
xy_dnp=nan(numel(X),numel(Y));
xz_dnp=nan(numel(X),numel(Z));

% XY plane scan, Z=0
for n=1:numel(X)

    % Parallel inner loop
    parfor k=1:numel(Y) %#ok<*PFBNS>
        
        % Localise arrays
        localpar=parameters; localint=inter;
        
        % Set proton position
        localint.coordinates{1,1}=[X(n) Y(k) 0.0];  
        
        % Spinach housekeeping
        spin_system=create(sys,localint);
        spin_system=basis(spin_system,bas);
        
        % Electron operators
        Ex=operator(spin_system,'Lx','E');
        Ez=operator(spin_system,'Lz','E');

        % Detection state
        Nz=state(spin_system,'Lz','1H');
        
        % Thermal equilibrium state
        H0=hamiltonian(assume(spin_system,'labframe'),'left');
        rho_eq=equilibrium(spin_system,H0);
        
        % Hamiltonian and relaxation superoperator
        H=hamiltonian(assume(spin_system,'esr'));
        R=relaxation(spin_system);
        
        % Add microwave terms to the Hamiltonian
        H=H+parameters.mw_pwr*Ex+parameters.mw_off*Ez;
        
        % Run the time evolution
        rho=evolution(spin_system,H+1i*R,[],rho_eq,...
                      localpar.t_evol,1,'final');
        
        % Get the DNP amplitude
        xy_dnp(n,k)=real(Nz'*rho)/real(Nz'*rho_eq);
    
    end

    % Update the plot
    subplot(1,2,1); imagesc(X,Y,xy_dnp);  
    ktitle('$Z=0$ plane'); axis square;
    kxlabel('Proton Y coordinate, Angstrom');
    kylabel('Proton X coordinate, Angstrom');
    cb=colorbar; cb.TickLabelInterpreter='latex';
    cb.Label.String='$^{1}$H DNP @ 20 ms';
    cb.Label.FontSize=12; cb.Limits=[-200 1];
    cb.Label.Interpreter='latex'; drawnow;

end

% XZ plane scan, Y = 0
for n=1:numel(X)

    % Parallel inner loop
    parfor k=1:numel(Z)
        
        % Localise arrays
        localpar=parameters; localint=inter;
        
        % Set proton position
        localint.coordinates{1,1}=[X(n) 0.0 Z(k)];
        
        % Spinach housekeeping
        spin_system=create(sys,localint);
        spin_system=basis(spin_system,bas);
        
        % Electron operators
        Ex=(operator(spin_system,'L-','E')+...
            operator(spin_system,'L+','E'))/2;
        Ez=operator(spin_system,'Lz','E');

        % Detection state
        Nz=state(spin_system,'Lz','1H');
        
        % Thermal equilibrium state
        H0=hamiltonian(assume(spin_system,'labframe'),'left');
        rho_eq=equilibrium(spin_system,H0);
        
        % Hamiltonian and relaxation superoperator
        H=hamiltonian(assume(spin_system,'esr'));
        R=relaxation(spin_system);
        
        % Add microwave terms to the Hamiltonian
        H=H+parameters.mw_pwr*Ex+parameters.mw_off*Ez;
        
        % Run the time evolution
        rho=evolution(spin_system,H+1i*R,[],rho_eq,...
                      localpar.t_evol,1,'final');
        
        % Get the DNP amplitude
        xz_dnp(n,k)=real(Nz'*rho)/real(Nz'*rho_eq);
    
    end

    % Update the plot
    subplot(1,2,2); imagesc(X,Z,xz_dnp);
    ktitle('$Y=0$ plane'); axis square; 
    kxlabel('Proton Z coordinate, Angstrom');
    kylabel('Proton X coordinate, Angstrom');
    cb=colorbar; cb.TickLabelInterpreter='latex';
    cb.Label.String='$^{1}$H DNP @ 20 ms';
    cb.Label.FontSize=12; cb.Limits=[-200 1];
    cb.Label.Interpreter='latex'; drawnow;

end

end

