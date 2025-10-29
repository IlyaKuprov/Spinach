% Cross-polarisation experiment in the doubly rotating frame. A single
% nitrogen-15 and a single proton. Spinning powder simulation starting
% from the thermal equilibrium using the grid-free version of the Fok-
% ker-Planck formalism.
%
% Calculation time: minutes with a Tesla A100 GPU,
%                   much longer otherwise.
%
% ilya.kuprov@weizmann.ac.il
% aanevzor@ncsu.edu

function cp_contact_mas_nh_gridfree()

% System specification
sys.magnet=9.394;
sys.isotopes={'15N','1H'};
          
% Interactions
inter.zeeman.scalar={0.00 0.00};
inter.coordinates={[0.00 0.00 0.00]
                   [0.00 0.00 1.05]};
inter.temperature=298;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% This needs a GPU
% sys.enable={'gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relevant operators
Nx=operator(spin_system,'Lx','15N');
Ny=operator(spin_system,'Ly','15N');
Hx=operator(spin_system,'Lx','1H');
Hy=operator(spin_system,'Ly','1H');

% MAS parameters
parameters.rate=10000;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.max_rank=42;
parameters.spins={'1H','15N'};
parameters.irr_powers=[5e4*ones(1,100);
                       4e4*ones(1,100)];
parameters.irr_opers={Hy Nx};
parameters.exc_opers={Hx Ny};
parameters.needs={'iso_eq'};
parameters.coil=state(spin_system,'Lx','15N');
parameters.time_steps=1e-5*ones(1,100);

% Simulation
fid=gridfree(spin_system,@cp_contact_hard,parameters,'nmr');

% Plot the answer
time_axis=[0 cumsum(parameters.time_steps)];
figure(); plot(time_axis,real(fid)); kgrid;
kylabel('$S_{\rm{X}}$ expectation value on $^{15}N$');  
kxlabel('time, seconds'); xlim tight;

end

