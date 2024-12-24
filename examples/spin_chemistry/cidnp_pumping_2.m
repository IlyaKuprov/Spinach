% A simulation of Figure 2A in IK's paper on chemically amplified 
% NOEs (https://doi.org/10.1016/j.jmr.2004.01.011).
%
% Calculation time: seconds.
%
% ilya.kuprov@weizmann.ac.il

function cidnp_pumping_2()

% Magnet field
sys.magnet=14.1;

% Isotopes
sys.isotopes={'1H','19F'};

% Chemical shifts
inter.zeeman.scalar={0.0, 0.0};

% Chemical shift anisotropies (DFT)
inter.zeeman.eigs{1}=[0 0 0];
inter.zeeman.euler{1}=[0 0 0];
inter.zeeman.eigs{2}=[-47 -16  63];
inter.zeeman.euler{2}=[0 0 0];

% Coordinates (DFT)
inter.coordinates={[0 0 0],[0 2.60 0]};

% J-coupling (expt)
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=50;

% Relaxation theories
inter.relaxation={'redfield'};
inter.rlx_keep='secular';
inter.tau_c={110e-12};
inter.equilibrium='IME';
inter.temperature=298;

% Formalism and basis
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get the Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Get relaxation matrix
R=relaxation(spin_system);

% Get detection states
Hz=state(spin_system,{'Lz'},{1}); 
Fz=state(spin_system,{'Lz'},{2}); 
HzFz=state(spin_system,{'Lz','Lz'},{1,2});
coil=[Fz -HzFz Hz];

% Add pumping terms
R=magpump(spin_system,R,2*Hz,1.3);
R=magpump(spin_system,R,2*Fz,34.0);

% Proton also has other relaxation
R=R-3.0*(Hz*Hz');

% Start at equilibrium
rho0=unit_state(spin_system)+2*Hz+2*Fz;

% Run the evolution
answer=evolution(spin_system,H+1i*R,coil,rho0,0.1,40,'multichannel');

% Do the plotting
figure(); scale_figure([1.5 1.0]);
time_axis=linspace(0,4,41);
subplot(1,3,1); plot(time_axis,real(answer(1,:))); 
ktitle('$H_{\mathrm{Z}}$'); kxlabel('time, s'); kgrid;
subplot(1,3,2); plot(time_axis,real(answer(2,:)));
ktitle('$F_{\mathrm{Z}}$'); kxlabel('time, s'); kgrid;
subplot(1,3,3); plot(time_axis,real(answer(3,:)));
ktitle('$H_{\mathrm{Z}}F_{\mathrm{Z}}$'); 
kxlabel('time, s'); kgrid;

end

