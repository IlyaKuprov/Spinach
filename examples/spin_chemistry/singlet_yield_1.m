% Liquid state magnetic field effect simulation on a radical
% pair with four nuclei using exponential recombination ki-
% netics model.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il
% h.j.hogben@chem.ox.ac.uk
% peter.hore@chem.ox.ac.uk

function singlet_yield_1()

% Unit magnet (field sweep)
sys.magnet=1;

% System specification
sys.isotopes ={'E','E','1H','1H','1H','1H'};
inter.zeeman.scalar={2.002 2.002 0 0 0 0};
inter.coupling.scalar = num2cell(mt2hz([0       0     0.195  0.195  0    0
                                        0       0     0      0     -1.3  0.2
                                        0.195   0     0      0      0    0  
                                        0.195   0     0      0      0    0  
                                        0      -1.3   0      0      0    0     
                                        0       0.2   0      0      0    0  ]/2));
                                    
% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=0;
bas.sym_spins={[3 4]};
bas.sym_group={'S2'};

% Fields and kinetics parameters
parameters.rates=[0.176 0.880 1.76 3.52 8.8 17.6 35.2 52.8]*1e6;
parameters.fields=1e-3*(0:0.01:5);
parameters.electrons=[1 2];
parameters.spins={'E'};
parameters.needs={'zeeman_op'};

% Disable ZTE
sys.disable={'zte'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
M=liquid(spin_system,@rydmr_exp,parameters,'labframe');

% Plot the answer
figure(); plot(parameters.fields,M,'r-'); 
kxlabel('magnetic field, Tesla'); kgrid;
kylabel('singlet recombination yield');

end

