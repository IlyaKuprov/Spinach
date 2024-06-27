% Singlet yield anisotropy calculation for a model radical pair
% reaction, Haberkorn recombination model.
%
% Run time: minutes on NVidia Titan V card, hours on CPU.
%
% i.kuprov@soton.ac.uk
% h.j.hogben@chem.ox.ac.uk
% peter.hore@chem.ox.ac.uk

function singlet_yield_anisotropy_3()

% Earth field
sys.magnet=50e-6;

% Isotopes
sys.isotopes={'E','E','14N','14N','1H','1H','1H'};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=5;

% Hyperfine coupling tensors
inter.coupling.matrix=cell(7);
inter.coupling.matrix{1,3}=mt2hz([-0.0989  0.0039  0.0000;  
                                   0.0039 -0.0989  0.0000;
                                   0.0000  0.0000  1.7569]);
inter.coupling.matrix{2,4}=mt2hz([-0.0336  0.0924 -0.1354;  
                                   0.0924  0.3303 -0.5318; 
                                  -0.1354 -0.5318  0.6680]);
inter.coupling.matrix{2,5}=mt2hz([-0.9920 -0.2091 -0.2003; 
                                  -0.2091 -0.2631  0.2803;
                                  -0.2003  0.2803 -0.5398]);
inter.coupling.matrix{2,6}=mt2hz([-0.9920 -0.2091 -0.2003;
                                  -0.2091 -0.2631  0.2803;
                                  -0.2003  0.2803 -0.5398]);
inter.coupling.matrix{2,7}=mt2hz([-0.9920 -0.2091 -0.2003;
                                  -0.2091 -0.2631  0.2803;
                                  -0.2003  0.2803 -0.5398]);

% Zeeman interactions
inter.zeeman.scalar={2.0023 2.0025 0 0 0 0 0};

% Kinetics parameters
inter.chem.rp_theory='haberkorn';
inter.chem.rp_electrons=[1 2];
inter.chem.rp_rates=[1e6 1e6];

% Sequence parameters
parameters.grid='leb_1ang_rank_63';
parameters.spins={'E'};
parameters.tol=1e-2;
parameters.verbose=0;
parameters.sum_up=0;

% Enable GPU arithmetic
sys.enable={'gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Run a simulation
[yields,grid]=powder(spin_system,@rydmr,parameters,'labframe');

% Do the plotting
figure(); plot(grid.betas,cell2mat(yields));
kxlabel('beta spherical angle, radians');
kylabel('singlet yield'); kgrid; axis tight;

end

