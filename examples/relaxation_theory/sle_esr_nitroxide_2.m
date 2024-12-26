% Slow motion regime simulation of an ESR spectrum of a nitroxide radical.
% Set to reproduce Figure 2 from the paper by Concilio et al.:
% 
%                     http://arxiv.org/abs/1511.01667
%
% Calculation time: seconds
%
% mariagrazia.concilio@sjtu.edu.cn
% ilya.kuprov@weizmann.ac.il

function sle_esr_nitroxide_2()

% Magnet field
sys.magnet=0.3343;

% Isotopes
sys.isotopes={'14N','E'};

% Coupling Matrices
inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=1e6*gauss2mhz([1.000e+001 3.544e+000 1.170e+001;
                                          3.544e+000 1.800e+001 5.072e+000;
                                          1.170e+001 5.072e+000 3.000e+001]);

% Zeeman Interactions
inter.zeeman.matrix=cell(1, 2);
inter.zeeman.matrix{2}=[2.0065794 -0.0007548 -0.0032848; 
                       -0.0007548  2.0056940 -0.0006008; 
                       -0.0032848 -0.0006008  2.0048920];

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% SLE parameters
parameters.max_rank=7;
parameters.tau_c=17e-9;
parameters.rho0=state(spin_system,'L+','E');
parameters.coil=state(spin_system,'L+','E');
parameters.decouple={};
parameters.spins={'E'};
parameters.sweep=[-2.2e8 2e8];
parameters.npoints=1650;
parameters.zerofill=1650;
parameters.axis_units='GHz-labframe';
parameters.invert_axis=1;
parameters.derivative=1;

% SLE simulation
spectrum_sle=gridfree(spin_system,@slowpass,parameters,'esr');

% SLE plotting
figure(); plot_1d(spin_system,real(spectrum_sle),parameters);

end

