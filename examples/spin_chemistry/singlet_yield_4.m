% Figure 3 from the paper by Till, Timmel, Brocklehurst and Hore:
%
%        http://dx.doi.org/10.1016/S0009-2614(98)01158-0
%
% Note: the original paper only uses electron Zeeman operators for
%       the field sweep, and therefore misses the effects associa-
%       ted with the rise in the nuclear Zeeman interaction on the
%       high field side of the resulting plot.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk
% h.j.hogben@chem.ox.ac.uk
% peter.hore@chem.ox.ac.uk

function singlet_yield_4()

% Unit magnet (field sweep)
sys.magnet=1;

% Spin system
sys.isotopes={'E','E','1H','1H'};

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Couplings
inter.zeeman.scalar={2.0023 2.0044 0 0};
inter.coupling.scalar{1,3}=gauss2mhz(35)*1e6;
inter.coupling.scalar{1,4}=gauss2mhz(30)*1e6;
inter.coupling.scalar{4,4}=0;

% Sequence parameters
parameters.fields=1e-3*10.^linspace(-5,3,2000); 
parameters.rates=[0.1 1.0 10.0 100.0 1000.0]*1e6;
parameters.electrons=[1 2];
parameters.spins={'E'};
parameters.needs={'zeeman_op'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
M=liquid(spin_system,@rydmr_exp,parameters,'labframe');

% Plotting
figure(); plot(linspace(-5,3,2000),M); kgrid;
kylabel('singlet recombination yield');
kxlabel('log(magnetic induction / mT)');

end

