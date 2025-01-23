% A scan through the microwave frequency range in a steady
% state DNP experiment for a single 15N labelled urea mole-
% cule at a specific orientation and a specific distance
% from a single electron.
%
% Laboratory frame DNP simulation is carried out with state
% space restriction to four-spin orders and a Weizmann DNP
% relaxation superoperator accounting for T1 and T2 and di-
% polar relaxation processes. Single crystal calculation.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il
% alexander.karabanov@nottingham.ac.uk
% walter.kockenberger@nottingham.ac.uk
% yonatan.hovav@weizmann.ac.il
% shimon.vega@weizmann.ac.il

function solid_effect_freq_scan_1()

% Magnetic field
sys.magnet=3.4;

% Spin system
sys.isotopes={'E','15N','1H','1H','15N','1H','1H'};
inter.coordinates={[ 0.00000000    0.00000000   10.14358975];
                   [-0.07640311    1.16112702   -0.61556225];
                   [ 0.08533754    1.99241453   -0.06489225];
                   [ 0.38824423    1.16155815   -1.51333625];
                   [ 0.07640311   -1.16112702   -0.61556225];
                   [-0.08533754   -1.99241453   -0.06489225];
                   [-0.38824423   -1.16155815   -1.51333625]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0'; 
bas.level=4;
bas.projections=[-2 -1 0 +1 +2];

% Relaxation theory
inter.relaxation={'weizmann'};
inter.rlx_keep='secular';
inter.equilibrium='zero';
inter.weiz_r1e=1e2;
inter.weiz_r1n=0.1;
inter.weiz_r2e=1e5;
inter.weiz_r2n=1e3;
inter.weiz_r1d=1e-3*ones(7,7);
inter.weiz_r2d=1e-3*ones(7,7);
inter.temperature=4.2;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E'};
parameters.mw_pwr=2*pi*100e3;
parameters.mw_frq=2*pi*[linspace(144.0,145.5,100)...
                        linspace(14.0,15.5,100)]*1e6;
parameters.coil=[state(spin_system,'Lz','1H')...
                 state(spin_system,'Lz','15N')];
parameters.mw_oper=(operator(spin_system,'L-','E')+...
                    operator(spin_system,'L+','E'))/2;
parameters.ez_oper=operator(spin_system,'Lz','E');
parameters.orientation=[pi/4 pi/5 pi/6];
parameters.method='lvn-backs';
parameters.needs={'aniso_eq'};
parameters.g_ref=spin_system.tols.freeg;

% Steady state simulation
answer=crystal(spin_system,@dnp_freq_scan,parameters,'esr');

% Plotting
figure(); scale_figure([1.75 1.75]);
subplot(2,2,1); plot(linspace(144.0,145.5,100),real(answer(1:100,1))); 
kxlabel('Microwave frequency, MHz'); kgrid; axis('tight');
kylabel('$S_\textrm{z}$ expectation value on $^{1}$H'); 
subplot(2,2,2); plot(linspace(14.0,15.5,100),real(answer(101:200,1))); 
kxlabel('Microwave frequency, MHz'); kgrid; axis('tight');
kylabel('$S_\textrm{z}$ expectation value on $^{1}$H'); 
subplot(2,2,3); plot(linspace(144.0,145.5,100),real(answer(1:100,2))); 
kxlabel('Microwave frequency, MHz'); kgrid; axis('tight');
kylabel('$S_\textrm{z}$ expectation value on $^{15}$N'); 
subplot(2,2,4); plot(linspace(14.0,15.5,100),real(answer(101:200,2))); 
kxlabel('Microwave frequency, MHz'); kgrid; axis('tight');
kylabel('$S_\textrm{z}$ expectation value on $^{15}$N'); 

end

