% A simple inversion-recovery experiment; longitudinal
% magnetisation is monitored as a function of time.
%
% Calculation time: seconds.
%
% ilya.kuprov@weizmann.ac.il

function inv_rec_1()

% Spin system
sys.magnet=14.1;
sys.isotopes={'1H'};

% Zeeman interactions
inter.zeeman.scalar={1.5};

% Complete basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
               
% Relaxation theory
inter.relaxation={'t1_t2'};
inter.r1_rates={5.0};
inter.r2_rates={5.0};
inter.equilibrium='dibari';
inter.rlx_keep='secular';
inter.temperature=298;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Isotropic thermal equilibrium
rho=equilibrium(spin_system);

% Detection state
coil=state(spin_system,'Lz','1H');

% Static Liouvillian superoperator
L=hamiltonian(assume(spin_system,'nmr'));

% Pulse operator
Lx=operator(spin_system,'Lx','1H');

% Relaxation superoperator
R=relaxation(spin_system);

% Inversion pulse
rho=step(spin_system,Lx,rho,pi);

% Evolution (1000 steps, 1.0 ms each)
answer=evolution(spin_system,L+1i*R,coil,rho,1e-3,1000,'observable');

% Plotting
kfigure(); scale_figure([1.00 0.65])
x_axis=linspace(0,1,1001); 
plot(x_axis,real(answer)); 
kylabel('$S_{\rm{Z}}$ expectation value'); 
kxlabel('time, seconds'); 
xlim tight; ylim padded; kgrid;

end

