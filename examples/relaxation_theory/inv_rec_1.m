% A simple inversion-recovery experiment.
%
% Calculation time: seconds.
%
% i.kuprov@soton.ac.uk

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

% Initial state - thermal equilibrium
[H,Q]=hamiltonian(assume(spin_system,'labframe'),'left');
rho=equilibrium(spin_system,H,Q,[0 0 0]);

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
figure(); x_axis=linspace(0,1,1001); 
plot(x_axis,real(answer)); 
kylabel('$S_{\rm{Z}}$ expectation value'); 
kxlabel('time, seconds'); 
xlim tight; kgrid;

end

