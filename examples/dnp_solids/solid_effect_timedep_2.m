% A simulation of solid effect DNP for a tilted linear chain of protons
% positioned at distances of 4+n^2 Angstrom with n=2:6 and the electron 
% located at the origin. Weizmann DNP relaxation model is used with se-
% cond order Krylov-Bogolyubov average Hamiltonian theory and state 
% space restriction to five-spin orders.
%
% Calculation time: minutes with a Tesla A100 GPU
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
% alexander.karabanov@nottingham.ac.uk
% walter.kockenberger@nottingham.ac.uk
% yonatan.hovav@weizmann.ac.il
% shimon.vega@weizmann.ac.il

function solid_effect_timedep_2()

% Magnetic field
sys.magnet=3.4;

% Electron
sys.isotopes{1}='E';
inter.coordinates{1}=[0 0 0];

% Nuclei
R=euler2dcm(pi/6,pi/7,pi/8);
for n=2:7
    sys.isotopes{n}='1H';
    inter.coordinates{n}=[0 0 4+n^2]*R;
end

% Relaxation theory
inter.relaxation={'weizmann'};
inter.weiz_r1e=1e2;
inter.weiz_r1n=0.1;
inter.weiz_r2e=1e5;
inter.weiz_r2n=1e3;
inter.weiz_r1d=zeros(7,7);
inter.weiz_r2d=zeros(7,7);
for n=2:6
    inter.weiz_r1d(n,n+1)=0.1;
    inter.weiz_r1d(n+1,n)=0.1;
    inter.weiz_r2d(n,n+1)=0.1;
    inter.weiz_r2d(n+1,n)=0.1;
end
inter.temperature=4.2;
inter.equilibrium='IME';
inter.rlx_keep='secular';

% Microwave power and offset
parameters.mw_pwr=2*pi*250e3;
parameters.nuclear_frq=2*pi*144.76e6;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=5;
bas.projections=[+2 +1 0 -1 -2];

% Algorithmic options
sys.disable={'krylov'};
sys.enable={'gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.theory='kb_second_order';
parameters.time_step=0.01;
parameters.n_steps=1000;

% Time domain simulation
parameters.calc_type='time_dependence';
answer=solid_effect(spin_system,parameters);

% Time domain plotting
time_axis=linspace(0,parameters.time_step*parameters.n_steps,...
                     parameters.n_steps+1);
figure(); scale_figure([1.75 0.75]);
subplot(1,2,1); plot(time_axis,real(answer(1,:)));
kylabel('$S_\textrm{z}$ expectation value'); 
kxlabel('time, seconds'); kgrid;
legend({'electron'},'Location','SouthEast'); 
set(gca,'XScale','log'); axis tight; 
subplot(1,2,2); plot(time_axis,real(answer(2:end,:))); 
kylabel('$S_\textrm{z}$ expectation value'); 
kxlabel('time, seconds'); kgrid;
set(gca,'XScale','log'); axis tight;
legend({'Proton 1','Proton 2','Proton 3',...
        'Proton 4','Proton 5','Proton 6'},...
        'Location','NorthWest'); drawnow;

% Steady state simulation
parameters.calc_type='steady_state';
answer=solid_effect(spin_system,parameters);
disp('Steady-state Tr(Sz*rho) on all spins:');
disp(real(answer));

end

