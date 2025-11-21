% Comparison between nitroxide simulation using SLE formalism
% and Redfield relaxation theory. 
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function sle_esr_nitroxide_1()

% Spin system properties
options.no_xyz=1;
[sys,inter]=g2spinach(gparse('../standard_systems/nitroxide.log'),...
                           {{'E','E'},{'N','14N'}},[0 0],options);
% Magnet induction
sys.magnet=3.5;

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
parameters.tau_c=5e-11;
parameters.rho0=state(spin_system,'L+','E');
parameters.coil=state(spin_system,'L+','E');
parameters.decouple={};
parameters.spins={'E'};
parameters.sweep=[-3e8 -1e8];
parameters.npoints=1024;
parameters.zerofill=1024;
parameters.axis_units='GHz-labframe';
parameters.invert_axis=1;

% SLE simulation
spectrum_sle=gridfree(spin_system,@slowpass,parameters,'esr');
spectrum_sle=fdvec(spectrum_sle,5,1);

% SLE plotting
kfigure(); subplot(1,2,1);
plot_1d(spin_system,real(spectrum_sle),parameters);
ktitle('SLE');

% BRW parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={5e-11};

% BRW housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% BRW parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'L+','E');
parameters.coil=state(spin_system,'L+','E');
parameters.sweep=[-3e8 -1e8];
parameters.npoints=1024;
parameters.zerofill=1024;
parameters.axis_units='GHz-labframe';

% BRW simulation
spectrum_brw=liquid(spin_system,@slowpass,parameters,'esr');
spectrum_brw=fdvec(spectrum_brw,5,1);

% BRW plotting
subplot(1,2,2);
plot_1d(spin_system,real(spectrum_brw),parameters);
ktitle('BRW');

end

