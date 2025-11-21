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

function tau_scan_si_sys_b()

% Magnetic field, Tesla
sys.magnet=14.1;

% Spin system
sys.isotopes={'1H','E','E'};
    
% Zeeman interactions
inter.zeeman.eigs{1,1}=[0 10 20];
inter.zeeman.euler{1,1}=[0 0 0];   
inter.zeeman.eigs{1,2}=[1.9778 1.9776 1.9776];          
inter.zeeman.euler{1,2}=[-0.180 0.017 0.194];     
inter.zeeman.eigs{1,3}=[2.0068 2.0038 2.0038];  
inter.zeeman.euler{1,3}=[0.632 0.783 1.086];

% Exchange coupling
inter.coupling.scalar=cell(3,3);
inter.coupling.scalar{2,3}=5e6;        

% coordinates
inter.coordinates={[ 0.000  0.000  0.000];
                   [ 6.000  0.030  0.317];
                   [-6.000 -0.038  0.535]};  

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

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
parameters.mw_frq=2*pi*linspace(-15,15,512)*1e6;

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

end

