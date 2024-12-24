% Overhauser type DNP in liquid phase at room temperature, using a continu-
% ous on-resonance CW irradiation of the electron ESR signal. The simulati-
% on uses Redfield theory to account for the dipolar cross-relaxation.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function odnp_liquid_1()

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
               
% Complete basis set
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

% Experiment paramaters
parameters.spins={'E'};
parameters.needs={'rho_eq'};
parameters.coil=[state(spin_system,{'Lz'},{1})...
                 state(spin_system,{'Lz'},{2})...
                 state(spin_system,{'Lz'},{3})];
parameters.mw_pwr=2*pi*1e6;
parameters.mw_off=0;
parameters.mw_oper=operator(spin_system,'Lx','E');
parameters.ez_oper=operator(spin_system,'Lz','E');
parameters.dt=1e-6;
parameters.nsteps=1e3;

% Simulation
answer=liquid(spin_system,@dnp_time_dep,parameters,'esr');

% Plotting
figure(); x_axis=linspace(0,1000,1001);
subplot(2,1,1); plot(x_axis,real(answer(3,:))); 
kxlabel('time, microseconds'); kgrid;
kylabel('$\langle E_{\rm{Z}} \rangle$');
subplot(2,1,2); plot(x_axis,real(answer(1:2,:))); 
kxlabel('time, microseconds'); kgrid;
kylabel('$\langle H_{\rm{Z}} \rangle$');
legend({'0.5 Angstrom from electron',...
        '1.5 Angstrom from electron'},...
        'Location','NorthEast');

end

