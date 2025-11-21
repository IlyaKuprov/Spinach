% Steady state populations of various spin states a function of rotational
% correlation time in a DNP experiment with two electrons connected by ex-
% change coupling, both coupled to a nucleus by dipolar couplings. Further
% particulars in:
%
%               https://doi.org/10.1016/j.jmr.2021.106940
%
% Calculation time: seconds
%
% mariagrazia.concilio@sjtu.edu.cn 
% ilya.kuprov@weizmann.ac.il

function states_vs_tauc_si_sys_a()

% Magnetic field, Tesla
sys.magnet=14.1;

% Spin system
sys.isotopes={'1H','E','E'};
      
% Zeeman interactions
inter.zeeman.eigs{1,1}=[0 10 20];
inter.zeeman.euler{1,1}=[0 0 0];
inter.zeeman.eigs{1,2}=[1.977873 1.977798 1.977792];     
inter.zeeman.euler{1,2}=[0 0 0];
inter.zeeman.eigs{1,3}=[1.977919 1.978000 1.978000];     
inter.zeeman.euler{1,3}=[-0.59 -0.10 0.49];

% Exchange coupling
inter.coupling.scalar=cell(3,3);
inter.coupling.scalar{2,3}=6.2e6;    

% Coordinates for anisotropic HF
inter.coordinates={[ 0.00  0.0000  0.0000];
                   [ 7.03  0.0187  0.9820];
                   [-7.03  0.2051 -1.0001]};   
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
parameters.mw_frq=2*pi*15.4e6;
parameters.g_ref=mean(inter.zeeman.eigs{1,2});  

% Correlation time grid
tau_c=linspace(50e-12,500e-12,64);

% Preallocate results
answer=zeros(12,numel(tau_c),'like',1i);

% Disable excessive printing and checking
sys.disable={'hygiene'}; sys.output='hush';

% Loop over correlation time
parfor n=1:numel(tau_c)

    % Localise parameter arrays
    locint=inter; locpar=parameters;
                            
    % Set correlation time      
    locint.tau_c={tau_c(n)};
                   
    % Spinach housekeeping
    spin_system=create(sys,locint);
    spin_system=basis(spin_system,bas);
    
   % Relevant states
    locpar.coil=[state(spin_system,{'L+'},{2})+2*state(spin_system,{'L+','Lz'},{2,3})...        
                 state(spin_system,{'L+'},{2})-2*state(spin_system,{'L+','Lz'},{2,3})...
                 state(spin_system,{'L+'},{3})+2*state(spin_system,{'Lz','L+'},{2,3})...        
                 state(spin_system,{'L+'},{3})-2*state(spin_system,{'Lz','L+'},{2,3})...
               2*state(spin_system,{'Lz','L+'},{1,2})...
               2*state(spin_system,{'Lz','L+'},{1,3})...
               4*state(spin_system,{'Lz','L+','Lz'},{1,2,3})...
               4*state(spin_system,{'Lz','Lz','L+'},{1,2,3})...
               2*state(spin_system,{'Lz','Lz'},{1,2})...
               2*state(spin_system,{'Lz','Lz'},{1,3})...
               4*state(spin_system,{'Lz','Lz','Lz'},{1,2,3})...
                 state(spin_system,{'Lz'},{1})];
    locpar.mw_oper=operator(spin_system,'Lx','E')/2;
    locpar.ez_oper=operator(spin_system,'Lz','E');

    % Steady state simulation
    answer(:,n)=liquid(spin_system,@dnp_freq_scan,locpar,'esr');
       
end

% Plotting
kfigure(); subplot(1,3,1); scale_figure([2.0 0.8]); 
plot(tau_c/1e-12',abs(answer(1,:))',...
     tau_c/1e-12',abs(answer(2,:))',...
     tau_c/1e-12',abs(answer(3,:))',...
     tau_c/1e-12',abs(answer(4,:))'); 
kxlabel('Correlation time, ps'); kgrid; 
kylabel('Steady state amplitude, a.u.');
legend({'${ \hat E_ {1+} + 2 \hat E_ {1+} \hat E_ {2z} }$',...
        '${ \hat E_ {1+} - 2 \hat E_ {1+} \hat E_ {2z} }$',...
        '${ \hat E_ {2+} + 2 \hat E_ {1z} \hat E_ {2+} }$',...
        '${ \hat E_ {2+} - 2 \hat E_ {1z} \hat E_ {2+} }$'},...
        'Interpreter','latex','Location','northeast');
xlim([min(tau_c) max(tau_c)]*1e12);
    
subplot(1,3,2); 
plot(tau_c/1e-12',abs(answer(5,:))',...
     tau_c/1e-12',abs(answer(6,:))',...
     tau_c/1e-12',abs(answer(7,:))',...
     tau_c/1e-12',abs(answer(8,:))'); 
kxlabel('Correlation time, ps'); kgrid; 
kylabel('Steady state amplitude, a.u.'); 
legend({'${ 2 \hat N_ {z} \hat E_ {1+} }$',...
        '${ 2 \hat N_ {z} \hat E_ {2+} }$',...
        '${ 4 \hat N_ {z} \hat E_ {1+} \hat E_ {2z} }$',...
        '${ 4 \hat N_ {z} \hat E_ {1z} \hat E_ {2+} }$'},...
        'Interpreter','latex','Location','northeast'); 
xlim([min(tau_c) max(tau_c)]*1e12);
    
subplot(1,3,3); 
plot(tau_c/1e-12',abs(answer(9 ,:))',...
     tau_c/1e-12',abs(answer(10,:))',...
     tau_c/1e-12',abs(answer(11,:))',...
     tau_c/1e-12',abs(answer(12,:))'); 
kxlabel('Correlation time, ps'); kgrid; 
kylabel('Steady state amplitude, a.u.'); 
legend({'${ 2 \hat N_ {z} \hat E_ {1z} }$',...
        '${ 2 \hat N_ {z} \hat E_ {2z} }$',...
        '${ 4 \hat N_ {z} \hat E_ {1z} \hat E_ {2z} }$',...
        '${ \hat N_ {z} }$'},...
        'Interpreter','latex','Location','northeast');
xlim([min(tau_c) max(tau_c)]*1e12);
   
end

