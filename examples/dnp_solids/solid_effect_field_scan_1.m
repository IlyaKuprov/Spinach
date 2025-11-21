% Magnetic field sweep DNP experiment involving a gadolinium ion,
% steady-state polarisation of a 15N nucleus is computed as a 
% function of the magnetic field offset.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il
% corzilius@em.uni-frankfurt.de

function solid_effect_field_scan_1()
   
% Magnetic field
sys.magnet=9.4509;

% Spin system
sys.isotopes={'E8','15N'};

% Electron g-tensor
inter.zeeman.eigs=cell(2,1);
inter.zeeman.euler=cell(2,1);
inter.zeeman.eigs{1}=[1.9918 1.9918 1.9918];
inter.zeeman.euler{1}=[0 0 0];

% Electron ZFS tensor
inter.coupling.eigs=cell(2,2);
inter.coupling.euler=cell(2,2);
inter.coupling.eigs{1,1}=570e6*[-1/3 -1/3 2*1/3];
inter.coupling.euler{1,1}=[0 0 0];

% Coordinates (Angstrom)
inter.coordinates={[0.00   0.00   0.00];
                   [3.00   0.00   0.00]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=[-2 -1 0 1 2];

% Relaxation theory
inter.relaxation={'t1_t2'};
inter.r1_rates={1e4 1e1};
inter.r2_rates={1e7 1e3};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.temperature=40.2;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E8'};
parameters.mw_pwr=1e5;
parameters.mw_frq=-14e8;
parameters.fields=linspace(-0.08,+0.08,512);
parameters.coil=state(spin_system,'Lz','15N');
parameters.mw_oper=(operator(spin_system,'L-','E8')+...
                    operator(spin_system,'L+','E8'))/2;
parameters.ez_oper=operator(spin_system,'Lz','E8');
parameters.grid='rep_2ang_6400pts_sph';
parameters.method='backslash';
parameters.needs={'aniso_eq'};

% Steady state simulation
spec=powder(spin_system,@dnp_field_scan,parameters,'esr');

% Plotting
kfigure(); plot(parameters.fields,real(spec)); kgrid;
axis tight; kxlabel('Magnetic field offset, Tesla');
kylabel('$S_\textrm{z}$ expectation value on $^{15}$N');

end

