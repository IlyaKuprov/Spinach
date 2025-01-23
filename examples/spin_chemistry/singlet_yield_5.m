% Figure 5 from the paper by Timmel, Till, Brocklehurst, McLauchlan
% and Hore:
%              http://dx.doi.org/10.1080/00268979809483134
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il
% h.j.hogben@chem.ox.ac.uk
% peter.hore@chem.ox.ac.uk

function singlet_yield_5()

% Unit magnet (field sweep)
sys.magnet=1;

% System specification
sys.isotopes={'E','E','1H','1H'};
inter.zeeman.scalar={2.0023 2.0023 0 0};
inter.coupling.scalar{1,3}=gauss2mhz(20)*1e6;
inter.coupling.scalar{2,4}=gauss2mhz(20)*1e6;
inter.coupling.scalar{4,4}=0; 

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Magnetic field and kinetics
parameters.fields=linspace(0,0.5*20/1e4,200);
parameters.rates=2*pi*[1e-4 1e-3 0.01 0.02 0.05 0.1 0.2]*gauss2mhz(20)*1e6;
parameters.electrons=[1 2];
parameters.spins={'E'};
parameters.needs={'zeeman_op'};

% Spinach run
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
M=liquid(spin_system,@rydmr_exp,parameters,'labframe');

% Plotting
figure(); plot(linspace(0,0.5,200),M); kgrid;
kylabel('singlet recombination yield');
kxlabel('$\omega/a$'); axis([0 0.5 0.46 0.58]);

end

