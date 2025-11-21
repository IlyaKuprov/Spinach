% Steady state nuclear magnetisation as a function of microwave frequency
% offset and the rotational correlation time in a DNP experiment with two
% electrons connected by exchange coupling, both coupled to a nucleus by
% dipolar couplings. Further particulars in:
%
%               https://doi.org/10.1016/j.jmr.2021.106940
%
% Calculation time: seconds
%
% mariagrazia.concilio@sjtu.edu.cn 
% ilya.kuprov@weizmann.ac.il

function tau_scan_main_text()

% Magnetic field, Tesla
sys.magnet=14.1;

% Spin system
sys.isotopes={'1H','E','E'};
      
% Zeeman interactions
inter.zeeman.eigs{1,1}=[0 10 20];
inter.zeeman.euler{1,1}=[0 0 0];   
inter.zeeman.eigs{1,2}=[2.0034 2.0038 2.0038];          
inter.zeeman.euler{1,2}=[-0.872 -0.013 0.868];     
inter.zeeman.eigs{1,3}=[2.0057 2.0030 2.0030];  
inter.zeeman.euler{1,3}=[-1.145 0.061 1.143];

% Exchange coupling
inter.coupling.scalar=cell(3,3);
inter.coupling.scalar{2,3}=3e6;        

% Cooridnates
inter.coordinates={[ 0.000  0.000  0.000];
                   [ 5.090  0.010  0.958];
                   [-5.090  0.061  1.032]};   
% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Turn off startup tests
sys.disable={'hygiene'};

% Relaxation theory
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.temperature=298;
sys.tols.rlx_integration=1e-10; % Needs to be this tight

% Sequence parameters
parameters.spins={'E'};
parameters.mw_pwr=2*pi*1e6;
parameters.method='lvn-backs';
parameters.needs={'rho_eq'};
parameters.g_ref=mean(inter.zeeman.eigs{1,2});

% Disable excessive printing and checking
sys.disable={'hygiene'}; sys.output='hush';

% Microwave frequency offset grid
parameters.mw_frq=2*pi*linspace(-5,10,512)*1e6;

% Correlation time grid
tau_c=linspace(50e-12,500e-12,128);

% Preallocate results
answer=zeros(numel(parameters.mw_frq),numel(tau_c),'like',1i);

% Loop over correlation time
parfor n=1:numel(tau_c)

    % Localise parameter arrays
    localpar=parameters; localint=inter;
                            
    % Set correlation time      
    localint.tau_c={tau_c(n)};
                   
    % Spinach housekeeping
    spin_system=create(sys,localint);
    spin_system=basis(spin_system,bas);
    
    % Relevant operators and states
    localpar.coil=state(spin_system,'Lz','1H');
    localpar.mw_oper=operator(spin_system,'Lx','E')/2;
    localpar.ez_oper=operator(spin_system,'Lz','E');

    % Steady state simulation
    answer(:,n)=liquid(spin_system,@dnp_freq_scan,localpar,'esr');

    % Thermal equilibium state
    H0=hamiltonian(assume(spin_system,'labframe'),'left');
    rho_eq=equilibrium(spin_system,H0);

    % Apply thermal equilibrium reference
    answer(:,n)=answer(:,n)/(localpar.coil'*rho_eq);
       
end

% Plotting
kfigure(); scale_figure([0.6 1.0]);
imagesc(tau_c/1e-12,(parameters.mw_frq/(2*pi*1e6)),real(answer));
kxlabel('$\tau_c$ / ps'); kcolourbar('Steady state $^{1}$H DNP');
kylabel('MW freq offset from $g_{iso}^{(1)}$, MHz');
text(80,-4,'A','FontSize',22,'Color','black','FontName','Arial');
   
end

