% Cavity-assisted DNP of a low-field Bi:Si donor, inspired by
% Fig. 1 and Appendix D of Lorriaux et al., arXiv:2610.19180v1.
% RF saturation prepares the lower donor manifold. Selective
% microwave and RF pulses redistribute its populations, then
% spin-cavity exchange and photon loss reset the donor.
%
% The full E/209Bi/Cn system uses sphten-liouv throughout.
% The cavity is disconnected during preparation and pulses.
% RF sweeps replace the experimental band-limited noise.
% Only the electron magnetic moment drives the pulses; the
% smaller nuclear magnetic moment is neglected. The complete
% isotropic hyperfine interaction is retained.
%
%                 result=cavity_dnp_minimal()
%
% result contains time in seconds, nuclear Iz, photon number,
% trace, and Fock-edge population as row vectors. It also holds
% preparation diagnostics, pulse frequencies (Hz) and durations
% (s), and the spin system, states, and operators.
%
% Calculation time: minutes; requests a 64-worker parallel pool
%
% ilya.kuprov@weizmann.ac.il

function result=cavity_dnp_minimal()

% Converged cavity cutoff, pulse rank, and saturation sweep count
nlevels=3;
max_rank=2;
ncycles=300;

% Magnet field and particles
sys.magnet=0.06537;
sys.isotopes={'E','209Bi',['C' num2str(nlevels)]};

% Donor interactions and temperature
inter.zeeman.scalar={2.00032,0,[]};
inter.coupling.scalar=cell(3,3);
inter.coupling.scalar{1,2}=1.475e9;
inter.modes.frqs={[] [] 0};
inter.temperature=0.025;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation={'none'};

% Algorithmic options
sys.output='hush';
sys.parallel={'local',64};
sys.enable={'paranoia'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system=assume(spin_system,'labframe');

% Thermal state and pulse operators
H=hamiltonian(spin_system);
H_left=hamiltonian(spin_system,'left');
rho=equilibrium(spin_system,H_left);
unit=state(spin_system,'E',1);
Lx=operator(spin_system,'Lx',1);
Ly=operator(spin_system,'Ly',1);

% Remove the disconnected cavity from the pulse calculation
blocks=reduce(spin_system,H+Lx+Ly,rho);
P=horzcat(blocks{:});
H=P'*H*P;
H_left=P'*H_left*P;
rho=P'*rho;
Lx=P'*Lx*P;
Ly=P'*Ly*P;
unit_donor=P'*unit;
Iz=P'*operator(spin_system,'Lz',2,'left')'*unit;

% Transition frequencies from the twenty donor energy levels
energies=sort(real(eig(full(H_left))));
levels=uniquetol(energies,size(H_left,1)*eps(norm(H_left,1)),DataScale=1);
assert(numel(levels)==20,'The full donor must have twenty distinct levels.');
rf_frqs=diff(levels(1:9))/(2*pi);
mw_frq=(levels(15)-levels(6))/(2*pi);
rf_frq=(levels(7)-levels(6))/(2*pi);

% Saturation pulse propagators for the eight lower-manifold lines
rf_maps=cell(8,1);
rf_amplitude=2*pi*100e3;
parfor line=1:8

    % Retain the full donor drift during each rectangular pulse
    [~,~,rf_maps{line}]=shaped_pulse_af(spin_system,H,Lx,0*Ly,rho,...
                                        rf_frqs(line),rf_amplitude,...
                                        10e-6,0,max_rank,'expm');

end

% Saturate the lower manifold with independent RF phase ensembles
for cycle=1:ncycles

    % Retain the preceding sweep for the saturation residual
    rho_previous=rho;
    for line=1:8

        % Average the donor phase after each saturation pulse
        rho=P'*coherence(spin_system,P*(rf_maps{line}*rho),{{[1 2],0}});

    end

end

% Record convergence and normalise the prepared donor state
scramble_change=norm(rho-rho_previous);
preparation_trace=real(unit_donor'*rho);
rho=rho/(unit_donor'*rho);
initial_iz=real(Iz'*rho);

% Selective microwave pulse followed by donor phase averaging
mw_time=17.02542402393068e-6;
rho_mw=shaped_pulse_af(spin_system,H,Lx,0*Ly,rho,mw_frq,...
                       rf_amplitude,mw_time,0,max_rank,'expm');
rho_mw=P'*coherence(spin_system,P*rho_mw,{{[1 2],0}});

% Selective RF pulse followed by donor phase averaging
rf_time=22.56805013451008e-6;
rho_rf=shaped_pulse_af(spin_system,H,Lx,0*Ly,rho_mw,rf_frq,...
                       rf_amplitude,rf_time,0,max_rank,'expm');
rho_rf=P'*coherence(spin_system,P*rho_rf,{{[1 2],0}});
pulse_iz=real(Iz'*[rho_mw rho_rf]);

% Attach cavity vacuum to the prepared and control states
vacuum=operator(spin_system,'BL1',3,'left');
rho_wait=vacuum*(P*rho_rf);
rho_no_rf=vacuum*(P*rho_mw);
rho_initial=vacuum*(P*rho);
rho_wait=rho_wait/(unit'*rho_wait);
rho_no_rf=rho_no_rf/(unit'*rho_no_rf);
rho_initial=rho_initial/(unit'*rho_initial);

% Connect the resonant cavity to a zero-temperature photon bath
inter.modes.frqs={[] [] mw_frq};
inter.modes.exchange=cell(3,3);
inter.modes.exchange{1,3}=100;
inter.modes.linewidths={[] [] 210e3};
inter.temperature=0;

% Rebuild the connected spin-cavity system
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system=assume(spin_system,'labframe');
spin_system.inter.modes.strength.exchange='rwa';
H=hamiltonian(spin_system);
R=relaxation(spin_system);
unit=state(spin_system,'E',1);

% Use Spinach's matrix-exponential branch for the reset propagator
spin_system.tols.small_matrix=size(H,1)+1;

% Nuclear polarisation, photon number, and cavity truncation monitor
Iz=operator(spin_system,'Lz',2,'left')'*unit;
number=operator(spin_system,'N',3,'left')'*unit;
fock_edge=operator(spin_system,['BL' num2str(nlevels)],3,'left')'*unit;
observables=[unit Iz number fock_edge];

% Cavity reset in the laboratory frame
time=0:0.001:20;
signals=evolution(spin_system,H+1i*R,observables,rho_wait,time(2),...
                  numel(time)-1,'multichannel');
signals=real(signals);

% Observable trajectories and net nuclear polarisation transfer
result.time=time;
result.trace=signals(1,:);
result.iz=signals(2,:);
result.photons=signals(3,:);
result.top_population=signals(4,:);
result.initial_iz=full(initial_iz);
result.pulse_iz=full(pulse_iz);
result.delta_iz=full(signals(2,end)-initial_iz);

% Preparation and pulse diagnostics
result.scramble_change=scramble_change;
result.preparation_trace=full(preparation_trace);
result.transition_hz=[mw_frq rf_frq];
result.saturation_hz=rf_frqs;
result.pulse_durations=[mw_time rf_time];

% Spin system and states underlying the trajectories
result.spin_system=spin_system;
result.H=H;
result.R=R;
result.observables=observables;
result.rho=rho_wait;
result.rho_no_rf=rho_no_rf;
result.rho_initial=rho_initial;
result.nlevels=nlevels;
result.max_rank=max_rank;
result.ncycles=ncycles;

end


