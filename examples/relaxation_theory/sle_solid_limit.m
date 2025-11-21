% Solid limit of Stochastic Liouville equation formalism.
%
% Calculation time: hours
%
% ilya.kuprov@weizmann.ac.il

function sle_solid_limit()

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
parameters.rho0=state(spin_system,'L+','E');
parameters.coil=state(spin_system,'L+','E');
parameters.decouple={};
parameters.spins={'E'};
parameters.sweep=[-2.2e8 2e8];
parameters.npoints=240;
parameters.zerofill=240;
parameters.axis_units='GHz-labframe';
parameters.invert_axis=1;
parameters.derivative=0;

% Ranks and correlation times
ranks=[3    7    15   30  ];
tau_c=[1e-9 1e-8 1e-7 1e-6];

% Start a figure
kfigure(); scale_figure([3.0 1.0]);

% Loop over correlation times
for n=1:numel(ranks)
    
    % Set the parameters
    parameters.max_rank=ranks(n);
    parameters.tau_c=tau_c(n);

    % SLE simulation
    spectrum_sle=gridfree(spin_system,@slowpass,parameters,'esr');

    % SLE plotting
    subplot(1,numel(ranks),n);
    plot_1d(spin_system,real(spectrum_sle),parameters);
    ktitle(['$\tau_{c}=10^{' num2str(log10(tau_c(n))) '}$']);
    kxlabel('El. Zeeman freq., GHz'); drawnow; 

end

end

