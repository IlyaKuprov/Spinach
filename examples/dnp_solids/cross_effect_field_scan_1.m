% Magnetic field sweep cross effect DNP experiment - steady-state proton 
% magnetisation under microwave iradiation as a function of the applied
% magnetic field. A powder average calculation.
%
% Calculation time: hours
%
% i.kuprov@soton.ac.uk

function cross_effect_field_scan_1()

% Magnetic field
sys.magnet=18.78;

% Spin system
sys.isotopes={'E','E','14N','1H'};

% Electron g-tensors
inter.zeeman.eigs=cell(4,1);
inter.zeeman.euler=cell(4,1);
inter.zeeman.eigs{1}=[2.0085 2.00605 2.00215];
inter.zeeman.euler{1}=[pi/2 pi/3 pi/4];
inter.zeeman.eigs{2}=[2.00319 2.00319 2.00258];
inter.zeeman.euler{2}=[pi/5 pi/6 pi/7];

% 14N quadrupolar tensor
inter.coupling.eigs=cell(4,4);
inter.coupling.euler=cell(4,4);
inter.coupling.eigs{3,3}=[-1e6 -1e6 2e6];
inter.coupling.euler{3,3}=[0 0 0];

% Coordinates (Angstrom)
inter.coordinates={[ 0.00   0.00   0.00];
                   [12.80   0.00   0.00];
                   []                      % 14N hyperfine specified below
                   [ 5.0    6.0    7.0]};
               
% Hyperfine couplings
inter.coupling.eigs{1,3}=[17.4e6 17.6e6 102e6];
inter.coupling.euler{1,3}=[0 0 0];

% Exchange coupling
inter.coupling.scalar=cell(4,4);
inter.coupling.scalar{1,2}=2*(-73e6);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'t1_t2'};
inter.r1_rates={1e5 1e5 1e4 1e3};
inter.r2_rates={1e7 1e7 1e5 1e4};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.temperature=10;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.mw_pwr=10e6;
parameters.mw_frq=0;
parameters.spins={'E'};
parameters.fields=linspace(-0.08,0.04,256);
parameters.coil=state(spin_system,'Lz','1H');
parameters.mw_oper=operator(spin_system,'Lx','E')/2;
parameters.ez_oper=operator(spin_system,'Lz','E');
parameters.grid='rep_2ang_1600pts_sph';
parameters.method='backslash';
parameters.needs={'aniso_eq'};

% Steady state simulation
answer=powder(spin_system,@dnp_field_scan,parameters,'esr');

% Plotting
figure(); plot(parameters.fields,real(answer)); kgrid;
axis tight; kxlabel('Magnetic field offset, Tesla');
kylabel('$S_\textrm{z}$ expectation value on $^{1}$H'); 


end

