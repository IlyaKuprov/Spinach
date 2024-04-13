% Cross-polarisation experiment in the doubly rotating frame. A single
% nitrogen-15 in a bath of 8 protons scattered on a 2 Angstrom radius
% sphere around it. Spinning powder simulation using a restricted Lio-
% uville space up to, and including, three-spin correlations.
%
% Calculation time: minutes on NVidia Tesla A100, much longer on CPU
%
% i.kuprov@soton.ac.uk
% aanevzor@ncsu.edu

function cp_contact_mas_nhh()

% System specification
sys.magnet=9.394;
sys.isotopes={'1H','1H','1H','1H',...
              '1H','1H','1H','1H','15N'};
          
% Interactions
inter.zeeman.scalar={0.6745  1.0368 -0.1495 -0.3171  ...
                     0.7233  0.4882  3.2662  0.0317  0.0000};
inter.coordinates={[-2.51887819   -0.99807636    0.87365165]
                   [-2.54376425    2.65945619   -1.20122894]
                   [-1.11551509    1.65289357   -1.19927242]
                   [-2.54058796    1.14568573   -2.07390040]
                   [-4.50217572    1.96708252    0.00023756]
                   [-4.50219487    0.45407982   -0.87377011]
                   [-2.54233875    1.14692134    2.07390117]
                   [-1.11551334    1.65103779    1.20034313]
                   [-2.67552180    0.95825426    0.00000000]};
inter.temperature=298;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=3;

% Algorithmic options
sys.tols.inter_cutoff=5.0;
sys.tols.prox_cutoff=4.0;
sys.disable={'trajlevel'};
sys.enable={'greedy','gpu'};

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
parameters.max_rank=3;
parameters.irr_powers=[5e4*ones(1,100);
                       4e4*ones(1,100)];
parameters.spins={'1H','15N'};
parameters.irr_opers={Hy Nx};
parameters.exc_opers={Hx Ny};
parameters.needs={'iso_eq'};
parameters.coil=state(spin_system,'Lx','15N');
parameters.time_steps=1e-5*ones(1,100);
parameters.grid='rep_2ang_100pts_sph';

% Simulation
fid=singlerot(spin_system,@cp_contact_hard,parameters,'nmr');

% Plot the answer
time_axis=[0 cumsum(parameters.time_steps)];
figure(); plot(time_axis,real(fid)); kgrid;
kylabel('$S_{\rm{X}}$ expectation value on $^{15}N$');  
kxlabel('time, seconds'); xlim tight;

end

