% Overhauser type DNP in liquid phase at room temperature, after a perfect
% inversion pulse on the electron ESR signal. The simulation uses Redfield
% theory to account for the dipolar cross-relaxation.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function odnp_liquid_2()

% Spin system
sys.magnet=3.4;
sys.isotopes={'1H','1H','E'};

% Zeeman interactions
inter.zeeman.matrix={[5 0 0; 0 5 0; 0 0 5]...
                     [5 0 0; 0 5 0; 0 0 5]...
                     [2.0023 0 0; 0 2.0025 0; 0 0 2.0027]};

% Coordinates (Angstrom)
inter.coordinates={[0.0 0.0 0.0]
                   [0.0 2.0 0.0]
                   [0.0 0.0 1.5]};
               
% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'redfield'};
inter.equilibrium='dibari';
inter.rlx_keep='secular';
inter.temperature=298;
inter.tau_c={10e-12};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Thermal equilibrium state
[H,Q]=hamiltonian(assume(spin_system,'labframe'),'left');
rho_eq=equilibrium(spin_system,H,Q,[0 0 0]);

% Electron control operator
Lx=operator(spin_system,'Lx','E');

% Inversion pulse on the electron
rho0=step(spin_system,Lx,rho_eq,pi);
             
% Experiment paramaters
parameters.spins={'E'};
parameters.rho0=rho0;
parameters.coil=[state(spin_system,{'Lz'},{1})...
                 state(spin_system,{'Lz'},{2})...
                 state(spin_system,{'Lz'},{3})];
parameters.mw_pwr=0;
parameters.mw_off=0;
parameters.mw_oper=operator(spin_system,'Lx','E');
parameters.ez_oper=operator(spin_system,'Lz','E');
parameters.dt=1e-6;
parameters.nsteps=1e3;

% Simulation
answer=liquid(spin_system,@dnp_time_dep,parameters,'esr');

% Plotting
kfigure(); x_axis=linspace(0,1000,1001);
subplot(2,1,1); plot(x_axis,real(answer(3,:))); 
kxlabel('time, microseconds'); kgrid;
kylabel('$\langle E_{\rm{Z}} \rangle$');
subplot(2,1,2); plot(x_axis,real(answer(1:2,:))); 
kxlabel('time, microseconds'); kgrid;
kylabel('$\langle H_{\rm{Z}} \rangle$');
legend({'0.5 Angstrom from electron',...
        '1.5 Angstrom from electron'},...
        'Location','SouthEast');

end

