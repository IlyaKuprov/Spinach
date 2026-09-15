% Efficiency of the z-filtered 27Al MQMAS pulse sequence with hard
% pulses and with the optimal control pulses produced by the other
% examples in this folder. Reproduces, using Spinach, the sequence
% efficiency calculation from
%
%           https://doi.org/10.26434/chemrxiv.15008427
%
% The nucleus has the quadrupolar coupling of aluminium acetylace-
% tonate (CQ=3.2 MHz, eta=0.16) and is spun at 12.5 kHz in a 400
% MHz magnet; the quadrupolar interaction is taken to second order
% in the rotating frame. The sequence is: excitation pulse, +MQ/-MQ
% coherence filter, conversion pulse, population filter, central
% transition selective pulse, and detection of the central transi-
% tion single-quantum coherence. The efficiency is the modulus of
% the detected element of the density matrix, normalised to the
% initial state as in the paper. The powder average runs over 400
% crystallite orientations at 32 initial rotor phases each. Hard
% pulses have the durations optimised in the paper. Optimal control
% waveforms are read from the files written by mq_excitation.m,
% mq_conversion.m, and ct_selective.m examples; the files supplied
% in this folder were produced by those examples on 128 cores.
%
% Calculation time: minutes on 128 cores.
%
% ilya.kuprov@weizmann.ac.il

function mqmas_efficiency()

% Coherence order, 3 or 5
mq_order=5;

% 400 MHz magnet
sys.magnet=2*pi*400e6/spin('1H');
sys.isotopes={'27Al'};

% Quadrupolar coupling and shielding anisotropy
inter.coupling.matrix{1,1}=eeqq2nqi(3.2e6,0.16,5/2,[0 0 0]);
inter.zeeman.eigs={[-5 -5 10]};
inter.zeeman.euler={[0 0 0]};

% Hilbert space formalism
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system=assume(spin_system,'labframe');

% Rotor phase resolved drift Hamiltonians for the whole sequence
parameters.spins={'27Al'};
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.grid='rep_2ang_400pts_sph';
parameters.n_ticks=160;
parameters.n_phases=32;
parameters.n_slices=1060;
drifts=mqmas_drifts(spin_system,parameters);
tick_dt=0.5e-6;

% Initial state, Iz
rho_init=state(spin_system,'Lz','27Al');
rho_init=rho_init/norm(rho_init,'fro');

% Control operators
Lx=operator(spin_system,'Lx','27Al');
Ly=operator(spin_system,'Ly','27Al');

% Hard pulse durations and amplitudes from the paper
if mq_order==3
    hard_durs=[4.2e-6 1.4e-6 9e-6];
else
    hard_durs=[4.4e-6 2.4e-6 9e-6];
end
hard_amps=2*pi*[100e3 100e3 9.3e3];

% Hard pulses sliced at the drift tick interval
hard_pulses=cell(1,3);
for k=1:3
    slice_durs=[tick_dt*ones(1,floor(hard_durs(k)/tick_dt)) mod(hard_durs(k),tick_dt)];
    slice_durs=slice_durs(slice_durs>1e-12);
    hard_pulses{k}={hard_amps(k)*[ones(1,numel(slice_durs)); zeros(1,numel(slice_durs))],slice_durs};
end

% Optimal control pulse sequence
exc=load(['mq_exc_' num2str(mq_order) 'q.mat'],'pulse','pulse_dt');
conv=load(['mq_conv_' num2str(mq_order) 'q.mat'],'pulse','pulse_dt');
ct=load('ct_pulse.mat','pulse','pulse_dt');
oc_pulses={{exc.pulse,exc.pulse_dt},{conv.pulse,conv.pulse_dt},{ct.pulse,ct.pulse_dt}};

% Loop over the two sequences
sequences={hard_pulses,oc_pulses}; efficiency=zeros(1,2);
for s=1:2

    % Parallel loop over the ensemble
    signals=zeros(1,numel(drifts)); pulses=sequences{s};
    parfor n=1:numel(drifts) %#ok<*PFBNS>

        % Start with the initial state at time zero
        rho=rho_init; current_time=0;

        % Loop over the three pulses
        for k=1:3

            % Loop over the slices of the pulse
            for m=1:numel(pulses{k}{2})

                % Drift Hamiltonian at the slice midpoint
                tick_idx=floor((current_time+pulses{k}{2}(m)/2)/tick_dt)+1;

                % Take a time step
                rho=step(spin_system,drifts{n}{tick_idx}+...
                         pulses{k}{1}(1,m)*Lx+pulses{k}{1}(2,m)*Ly,...
                         rho,pulses{k}{2}(m));
                current_time=current_time+pulses{k}{2}(m);

            end

            % Coherence filters after excitation and conversion
            if k==1
                rho=coherence(spin_system,rho,{{'27Al',[+mq_order -mq_order]}});
            elseif k==2
                rho=coherence(spin_system,rho,{{'27Al',0}});
            end

        end

        % Central transition single-quantum coherence
        signals(n)=rho(3,4);

    end

    % Powder average
    efficiency(s)=abs(mean(signals));

end

% Report the outcome
disp(['Hard pulse ' num2str(mq_order) 'QMAS efficiency:       ' num2str(efficiency(1))]);
disp(['Optimal control ' num2str(mq_order) 'QMAS efficiency: ' num2str(efficiency(2))]);
disp(['Signal enhancement factor:            ' num2str(efficiency(2)/efficiency(1))]);

end

