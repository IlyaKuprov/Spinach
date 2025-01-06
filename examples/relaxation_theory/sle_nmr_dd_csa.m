% 15N-1H DD-CSA cross-correlation in a protein amide bond 
% spin system using SLE formalism and Redfield relaxation
% theory.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function sle_nmr_dd_csa()

% System specification
sys.magnet=14.1;
sys.isotopes={'15N','1H'};
inter.zeeman.matrix{1}=[14.89    34.35     0.00
                        34.35   145.74     0.00
                         0.00     0.00   150.51];
inter.zeeman.matrix{2}=[30.75    1.31    0.00
                         1.31   22.65    0.00
                         0.00    0.00   11.80];
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=10;
inter.coordinates={[-0.451455   -0.678015    0.000000]
                   [-1.475290   -0.641823    0.000000]};

% Proximity cut-off
sys.tols.prox_cutoff=4.0;
               
% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% SLE housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% SLE parameters
parameters.max_rank=10;
parameters.tau_c=5e-9;
parameters.rho0=state(spin_system,'L+','15N');
parameters.coil=state(spin_system,'L+','15N');
parameters.decouple={};
parameters.spins={'15N'};
parameters.sweep=[-0.633e4 -0.629e4];
parameters.npoints=2048;
parameters.axis_units='Hz';

% SLE simulation
spectrum_sle=gridfree(spin_system,@slowpass,parameters,'nmr');

% SLE plotting
figure(); subplot(1,2,1); 
plot_1d(spin_system,real(spectrum_sle),parameters);
ktitle('SLE'); kylabel('amplitude, a.u.');

% BRW parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={5e-9};

% BRW housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% BRW parameters
parameters.spins={'15N'};
parameters.rho0=state(spin_system,'L+','15N');
parameters.coil=state(spin_system,'L+','15N');
parameters.sweep=[-0.633e4 -0.629e4];
parameters.npoints=2048;
parameters.axis_units='Hz';

% BRW simulation
spectrum_brw=liquid(spin_system,@slowpass,parameters,'nmr');

% BRW plotting
subplot(1,2,2); 
plot_1d(spin_system,real(spectrum_brw),parameters);
ktitle('BRW');

end

