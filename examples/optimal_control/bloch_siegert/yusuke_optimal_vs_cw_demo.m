% Bloch-Siegert-aware phase optimisation compared to a simple constant-
% phase low-power cycle. This is the control-side companion to
% yusuke_14n_broadening_demo.m: the task is a reduced-model surrogate for
% low-power offset-tolerant decoupling, formulated as an identity cycle
% that should preserve magnetisation across offset and B1 distributions.
%
% The "non-optimal pulse" is a constant-phase X pulse with the same RF
% amplitude and total duration. It behaves like a simple CW-style cycle:
% acceptable near the design point, but poor once offset and B1 errors
% are included. The "optimal pulse" is a phase-modulated waveform
% optimised with Bloch-Siegert corrections enabled.
%
% aditya.dev@weizmann.ac.il

function yusuke_optimal_vs_cw_demo()

rng(1);

% Magnetic field corresponding to 800 MHz 1H
sys.magnet=18.8;

% Single-spin surrogate model
sys.isotopes={'13C'};
inter.zeeman.scalar={0};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% States representing the identity map
Sx=state(spin_system,'Lx','13C'); Sx=Sx/norm(full(Sx),2);
Sy=state(spin_system,'Ly','13C'); Sy=Sy/norm(full(Sy),2);
Sz=state(spin_system,'Lz','13C'); Sz=Sz/norm(full(Sz),2);

% Control operators and drift
Lx=operator(spin_system,'Lx','13C');
Ly=operator(spin_system,'Ly','13C');
Lz=operator(spin_system,'Lz','13C');
H=hamiltonian(assume(spin_system,'nmr'));

% Practical low-power regime inspired by the 14N decoupling papers:
% 70 kHz MAS, tp ~ 10 us, and nu_14N in the 15-23 kHz range. In this
% reduced model, we use a 20 kHz carrier and a 10-step cycle of 10 us
% elements, giving a 4pi constant-phase reference pulse. The effective
% offset window below is intentionally narrower than the full laboratory-
% frame 14N dispersion because this script is a reduced control surrogate,
% not a full QJF/MAS quadrupolar simulation.
rf_nominal_hz=20e3;
nsteps=10;
pulse_dt=10e-6*ones(1,nsteps);
train_offsets_hz=linspace(-12e3,12e3,7);
train_b1_scales=[0.95 1.00 1.05];
eval_offsets_hz=linspace(-20e3,20e3,61);
eval_b1_scales=linspace(0.90,1.10,9);

% Define the phase-only optimisation problem
control.isotopes={'13C'};
control.channels=[1; 1];
control.drifts={{H}};
control.operators={Lx,Ly};
control.rho_init={Sx Sy Sz};
control.rho_targ={Sx Sy Sz};
control.off_ops={Lz};
control.offsets={train_offsets_hz};
control.pwr_levels=2*pi*rf_nominal_hz*train_b1_scales;
control.amplitudes=ones(1,nsteps);
control.pulse_dt=pulse_dt;
control.method='lbfgs';
control.max_iter=40;
control.plotting={};

% Enable Bloch-Siegert corrections in the optimiser and simulator
control.bsiegert=true();

% Spinach housekeeping
spin_system_bs=optimcon(spin_system,control);

% Optimise a phase-modulated identity cycle
guess=2*pi*rand(1,nsteps);
phi_opt=fmaxnewton(spin_system_bs,@grape_phase,guess);

% Constant-phase baseline with the same duration and RF power
phi_cw=zeros(1,nsteps);

% Evaluate both waveforms on a dense offset and B1 grid
[profile_opt,mean_offset_opt]=offset_profile(spin_system_bs,H,{Lx,Ly},...
    control.pulse_dt,Sx,Sy,Sz,Lz,phi_opt,rf_nominal_hz,eval_offsets_hz,...
    eval_b1_scales);

[profile_cw,mean_offset_cw]=offset_profile(spin_system_bs,H,{Lx,Ly},...
    control.pulse_dt,Sx,Sy,Sz,Lz,phi_cw,rf_nominal_hz,eval_offsets_hz,...
    eval_b1_scales);

[~,mean_b1_opt]=b1_profile(spin_system_bs,H,{Lx,Ly},...
    control.pulse_dt,Sx,Sy,Sz,Lz,phi_opt,rf_nominal_hz,eval_offsets_hz,...
    eval_b1_scales);

[~,mean_b1_cw]=b1_profile(spin_system_bs,H,{Lx,Ly},...
    control.pulse_dt,Sx,Sy,Sz,Lz,phi_cw,rf_nominal_hz,eval_offsets_hz,...
    eval_b1_scales);

% Training-set figures of merit
train_opt=ensemble_score(spin_system_bs,H,{Lx,Ly},control.pulse_dt,...
    Sx,Sy,Sz,Lz,phi_opt,rf_nominal_hz,train_offsets_hz,train_b1_scales);

train_cw=ensemble_score(spin_system_bs,H,{Lx,Ly},control.pulse_dt,...
    Sx,Sy,Sz,Lz,phi_cw,rf_nominal_hz,train_offsets_hz,train_b1_scales);

% Plot the result
figure;
time_axis_us=cumsum(control.pulse_dt)*1e6;

subplot(3,1,1);
plot(time_axis_us,unwrap(phi_opt),'LineWidth',1.2); hold on;
plot(time_axis_us,phi_cw,'LineWidth',1.2);
grid on;
xlabel('Time / \mus');
ylabel('Phase / rad');
title('Optimised low-power phase cycle versus constant-phase baseline');
legend({'Optimised','Constant-phase'},'Location','Best');

subplot(3,1,2);
plot(eval_offsets_hz/1e3,mean_offset_opt,'LineWidth',1.5); hold on;
plot(eval_offsets_hz/1e3,mean_offset_cw,'LineWidth',1.5);
grid on;
xlabel('Offset / kHz');
ylabel('Mean identity fidelity');
title('Average over B_1 distribution');
legend({'Optimised','Constant-phase'},'Location','SouthWest');

subplot(3,1,3);
plot(eval_b1_scales,mean_b1_opt,'LineWidth',1.5); hold on;
plot(eval_b1_scales,mean_b1_cw,'LineWidth',1.5);
grid on;
xlabel('B_1 scaling factor');
ylabel('Mean identity fidelity');
title('Average over offset distribution');
legend({'Optimised','Constant-phase'},'Location','SouthWest');

% Secondary robustness figures
figure;

subplot(2,1,1);
imagesc(eval_offsets_hz/1e3,eval_b1_scales,profile_opt);
axis xy;
colorbar;
xlabel('Offset / kHz');
ylabel('B_1 scaling factor');
title('Optimised pulse fidelity map');

subplot(2,1,2);
imagesc(eval_offsets_hz/1e3,eval_b1_scales,profile_cw);
axis xy;
colorbar;
xlabel('Offset / kHz');
ylabel('B_1 scaling factor');
title('Constant-phase pulse fidelity map');

% Console summary
fprintf('\n');
fprintf('Yusuke optimal-vs-CW control demo\n');
fprintf('  Training ensemble mean fidelity, optimised:      %.4f\n',mean(train_opt(:)));
fprintf('  Training ensemble mean fidelity, constant phase: %.4f\n',mean(train_cw(:)));
fprintf('  Training ensemble worst-case fidelity, optimised:      %.4f\n',min(train_opt(:)));
fprintf('  Training ensemble worst-case fidelity, constant phase: %.4f\n',min(train_cw(:)));
fprintf('\n');
fprintf(['Interpretation: the optimised low-power phase cycle is a surrogate ' ...
         'for offset-tolerant decoupling, while the constant-phase cycle ' ...
         'is the non-optimal baseline.\n']);

end

function [score_map,mean_profile]=offset_profile(spin_system,H,controls,...
    pulse_dt,Sx,Sy,Sz,Lz,phi_profile,rf_nominal_hz,offsets_hz,b1_scales)

% Fidelity map across offsets and B1 scales
score_map=zeros(numel(b1_scales),numel(offsets_hz));
for n=1:numel(offsets_hz)
    for k=1:numel(b1_scales)
        score_map(k,n)=single_score(spin_system,H,controls,pulse_dt,...
            Sx,Sy,Sz,Lz,phi_profile,rf_nominal_hz,offsets_hz(n),b1_scales(k));
    end
end
mean_profile=mean(score_map,1);

end

function [score_map,mean_profile]=b1_profile(spin_system,H,controls,...
    pulse_dt,Sx,Sy,Sz,Lz,phi_profile,rf_nominal_hz,offsets_hz,b1_scales)

% Fidelity map across offsets and B1 scales
score_map=zeros(numel(b1_scales),numel(offsets_hz));
for n=1:numel(offsets_hz)
    for k=1:numel(b1_scales)
        score_map(k,n)=single_score(spin_system,H,controls,pulse_dt,...
            Sx,Sy,Sz,Lz,phi_profile,rf_nominal_hz,offsets_hz(n),b1_scales(k));
    end
end
mean_profile=mean(score_map,2);

end

function score=ensemble_score(spin_system,H,controls,pulse_dt,Sx,Sy,Sz,...
    Lz,phi_profile,rf_nominal_hz,offsets_hz,b1_scales)

% Training-ensemble fidelity matrix
score=zeros(numel(b1_scales),numel(offsets_hz));
for n=1:numel(offsets_hz)
    for k=1:numel(b1_scales)
        score(k,n)=single_score(spin_system,H,controls,pulse_dt,...
            Sx,Sy,Sz,Lz,phi_profile,rf_nominal_hz,offsets_hz(n),b1_scales(k));
    end
end

end

function score=single_score(spin_system,H,controls,pulse_dt,Sx,Sy,Sz,...
    Lz,phi_profile,rf_nominal_hz,offset_hz,b1_scale)

% Build the physical waveform
amp_profile=2*pi*rf_nominal_hz*b1_scale*ones(size(phi_profile));
[CLx,CLy]=polar2cartesian(amp_profile,phi_profile);


% Apply the Bloch-Siegert-aware pulse
[controls_aug,pulse_aug]=bloch_siegert(spin_system,controls,...
                                             {CLx,CLy});
rho=shaped_pulse_xy(spin_system,H+2*pi*offset_hz*Lz,controls_aug,...
                    pulse_aug,pulse_dt,[Sx,Sy,Sz],'expv-pwc');

% Average identity fidelity across Cartesian basis states
score=real(trace([Sx, Sy, Sz]'*rho))/3;

end


