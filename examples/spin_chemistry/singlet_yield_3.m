% Figure 1 from the paper by Timmel, Till, Brocklehurst, McLauchlan
% and Hore:
%                http://dx.doi.org/10.1080/00268979809483134
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il
% h.j.hogben@chem.ox.ac.uk
% peter.hore@chem.ox.ac.uk

function singlet_yield_3()

% Unit magnet (field sweep)
sys.magnet=1;

% Spin system
sys.isotopes={'E','E','1H'};

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Couplings
inter.zeeman.scalar={2.0023 2.0023 0};
inter.coupling.scalar{1,3}=gauss2mhz(20)*1e6;
inter.coupling.scalar{3,3}=0;

% Sequence parameters
parameters.fields=linspace(0,3*20/1e4,200);
parameters.rates=2*pi*[0.005 0.02 0.05 0.1...
                      0.15 0.2 0.3 0.5 2.0]*gauss2mhz(20)*1e6;
parameters.electrons=[1 2];
parameters.spins={'E'};
parameters.needs={'zeeman_op'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
M=liquid(spin_system,@rydmr_exp,parameters,'labframe');

% Plotting
kfigure(); plot(linspace(0,3,200),M); kgrid;
kylabel('singlet recombination yield');
kxlabel('$\omega/a$'); axis([0 3 0.2 1.0]);

end

