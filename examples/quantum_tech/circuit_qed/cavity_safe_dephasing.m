% Stark-assisted flux-noise evasion (SAFE) for a cavity dispersively
% coupled to a flux-tunable transmon, Chapter 4 of Yunwei Lu's PhD
% thesis (Northwestern University, 2026). The transmon frequency ri-
% des on 1/f flux noise, which the dispersive coupling passes on to
% the cavity frequency; a weak off-resonant drive on the transmon
% Stark-shifts the cavity transition by an amount whose flux depen-
% dence cancels the direct one at the dynamical sweet spot. Flux
% noise trajectories are synthesised with pink_noise.m, each one is
% propagated under the Lindblad master equation in the frame rota-
% ting with the drive, and the ensemble average of the cavity cohe-
% rence between the dressed |g,0> and |g,1> states yields the pure
% dephasing time with and without the drive. The sweet spot drive
% amplitude is located numerically from the flux susceptibility of
% the dressed cavity transition and compared with Eq. (4.22) of the
% thesis; the noise spectrum and the coherence decay reproduce Figs.
% E.1 and 4.3(d) of the thesis with representative parameters.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function cavity_safe_dephasing()

% Transmon anharmonicity, transmon-cavity exchange coupling, and their detuning, Hz
anharm=-200e6; g_bc=100e6; delta_bc=1.0e9;

% Dispersive shift and the sensitivities of the cavity terms to the transmon frequency
chi=2*(g_bc/delta_bc)^2*anharm; sens_c=(g_bc/delta_bc)^2; sens_x=-4*anharm*g_bc^2/delta_bc^3;

% Transmon frequency sensitivity to flux, rad/s per flux quantum
dwb_dphi=2*pi*6e9;

% Flux noise amplitude in flux quanta, infrared and ultraviolet cutoffs, Hz
noise_amp=1e-5; f_ir=200; f_uv=5e7;

% Transmon-drive detuning, Hz
delta_bd=-19e6;

% Transmon and cavity relaxation times, seconds
t1_b=50e-6; t1_c=2.1e-3;

% Time step, number of steps, coherence sampling stride, and trajectory count
dt=1e-8; nsteps=25000; stride=100; ntraj=200;

% Length of the noise synthesis grid, long enough to resolve the infrared cutoff
nlong=2^19;

% Four-level transmon in the drive frame and a three-level cavity in its own frame
sys.magnet=0; sys.isotopes={'T4','C3'};
inter.modes.frqs={delta_bd 0};
inter.modes.anharms={anharm []};
inter.modes.lifetimes={t1_b t1_c};
inter.temperature=0;
bas.approximation='none';

% Hilbert space twin for the dressed states
bas.formalism='zeeman-hilb';
ss_hilb=create(sys,inter);
ss_hilb=basis(ss_hilb,bas);

% Hilbert space Hamiltonian, drive operator, and flux noise operator, Eqs. (4.16) and (4.33)
H_hilb=hamiltonian(assume(ss_hilb,'labframe'))+2*pi*chi*operator(ss_hilb,{'N','N'},{1,2});
D_hilb=operator(ss_hilb,'C',1)+operator(ss_hilb,'A',1);
N_hilb=operator(ss_hilb,'N',1)+sens_c*operator(ss_hilb,'N',2)+...
       sens_x*operator(ss_hilb,{'N','N'},{1,2});

% Positions of the bare |g,0> and |g,1> states
idx0=find(diag(state(ss_hilb,{'BL1','BL1'},{1,2}))>0.5);
idx1=find(diag(state(ss_hilb,{'BL1','BL2'},{1,2}))>0.5);

% Analytical sweet spot drive amplitude, Eq. (4.22), rad/s
omega0_an=2*pi*sqrt(delta_bd^3/(4*anharm));

% Numerical sweet spot from the flux susceptibility of the dressed cavity transition
omega0=fzero(@(omega)suscept(H_hilb,N_hilb,D_hilb,omega,idx0,idx1),[0.5 1.5]*omega0_an);
disp(['sweet spot drive amplitude, MHz: analytical ' num2str(omega0_an/(2*pi*1e6),'%.3f') ...
      ', numerical ' num2str(omega0/(2*pi*1e6),'%.3f')]);

% Liouville space for the dissipative dynamics
bas.formalism='zeeman-liouv';
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Drift, drive, and flux noise generators with the mode dissipators
L_drift=hamiltonian(assume(spin_system,'labframe'))+...
        2*pi*chi*operator(spin_system,{'N','N'},{1,2})+1i*relaxation(spin_system);
L_drive=operator(spin_system,'C',1)+operator(spin_system,'A',1);
L_noise=operator(spin_system,'N',1)+sens_c*operator(spin_system,'N',2)+...
        sens_x*operator(spin_system,{'N','N'},{1,2});

% Flux noise trajectories windowed from the long grid, transmon frequency units, rad/s
rng(1); dwb=zeros(ntraj,nsteps);
for m=1:ntraj
    dphi=pink_noise(noise_amp,dt,nlong,1,f_ir,f_uv); dwb(m,:)=dwb_dphi*dphi(1:nsteps);
end

% Frequency grid for the propagator table
ngrid=401; dw_max=max(abs(dwb(:))); dw_grid=linspace(-dw_max,dw_max,ngrid);
grid_idx=round((dwb+dw_max)*(ngrid-1)/(2*dw_max))+1;

% Time axis of the recorded coherence
time_axis=dt*stride*(1:(nsteps/stride));

% Preallocate the coherence decays and dephasing times
signals=zeros(2,numel(time_axis)); t2_times=zeros(1,2); drives=[0 omega0];

% Loop over the undriven and the driven case
for k=1:2

    % Dressed states adiabatically connected to |g,0> and |g,1>
    [V,E]=eig(full(H_hilb+drives(k)*D_hilb)); [~,k0]=max(abs(V(idx0,:))); [~,k1]=max(abs(V(idx1,:)));
    phi0=V(:,k0); phi1=V(:,k1); disp(['dressed transition frequency, kHz: ' num2str((E(k1,k1)-E(k0,k0))/(2*pi*1e3),'%.2f')]);

    % Equal superposition initial state and the coherence detection state
    rho=hilb2liouv((phi0+phi1)*(phi0+phi1)'/2,'statevec');
    coil=hilb2liouv(phi0*phi1','statevec');

    % Table of propagators over the frequency grid
    P=cell(1,ngrid);
    for n=1:ngrid
        P{n}=full(propagator(spin_system,L_drift+drives(k)*L_drive+dw_grid(n)*L_noise,dt));
    end

    % Ensemble average of the coherence over the noise trajectories
    coherence=zeros(1,numel(time_axis));
    for m=1:ntraj
        rho_traj=rho;
        for n=1:nsteps
            rho_traj=P{grid_idx(m,n)}*rho_traj;
            if mod(n,stride)==0
                coherence(n/stride)=coherence(n/stride)+coil'*rho_traj;
            end
        end
    end
    signals(k,:)=2*abs(coherence)/ntraj;

    % Coherence time from the 1/e crossing of the Gaussian decay or from an exponential fit
    if signals(k,end)<exp(-1)
        t2_times(k)=time_axis(find(signals(k,:)<exp(-1),1));
    else
        decay_fit=polyfit(time_axis,log(signals(k,:)),1); t2_times(k)=-1/decay_fit(1);
    end

end

% Pure dephasing times, Eq. (E.2), with the cavity photon loss removed
tphi_times=1./(1./t2_times-1/(2*t1_c));

% Analytical undriven pure dephasing time, Eq. (4.89), and the sweet spot residual, Eq. (4.87)
tphi_undriven=1/(noise_amp*dwb_dphi*sens_c*sqrt(2*abs(log(2*pi*f_ir*t2_times(1)))));
tphi_driven=1/((delta_bd/(4*anharm))*(dwb_dphi*noise_amp)^2/abs(delta_bd)+(delta_bd/(4*anharm))^2/t1_b);
disp(['pure dephasing time without the drive, us: simulated ' num2str(1e6*tphi_times(1),'%.1f') ...
      ', analytical ' num2str(1e6*tphi_undriven,'%.1f')]);
disp(['pure dephasing time with the drive, us: simulated ' num2str(1e6*tphi_times(2),'%.1f') ...
      ', analytical ' num2str(1e6*tphi_driven,'%.1f')]);

% Power spectral density of the windowed flux noise trajectories, Fig. E.1
freqs=(1:(nsteps/2))/(nsteps*dt); psd=mean(abs(fft(dwb/dwb_dphi,[],2)).^2,1)*dt/nsteps;
psd=psd(2:(nsteps/2+1)); band=(freqs>2e4)&(freqs<1e7);

% Validate the noise spectrum, the sweet spot location, and the dephasing suppression
if abs(mean(psd(band).*freqs(band))/noise_amp^2-1)>0.2
    error('flux noise spectrum deviates from 1/f.');
end
if abs(omega0/omega0_an-1)>0.5
    error('numerical sweet spot deviates from the analytical estimate.');
end
if (tphi_times(1)>2*tphi_undriven)||(tphi_times(1)<tphi_undriven/2)
    error('undriven dephasing time deviates from the analytical estimate.');
end
if tphi_times(2)<5*tphi_times(1)
    error('SAFE drive did not suppress the cavity dephasing.');
end

% Plot the noise spectrum and the coherence decays
kfigure(); subplot(1,2,1); loglog(freqs,psd,'.',freqs,noise_amp^2./freqs,'r-');
kgrid; kxlabel('frequency, Hz'); kylabel('$S(f)$, $\Phi_0^2$/Hz'); ktitle('flux noise spectrum');
subplot(1,2,2); plot(1e6*time_axis,signals'); kgrid; ylim([0 1.05]);
kxlabel('time, $\mu$s'); kylabel('cavity coherence'); ktitle('SAFE Ramsey decay');
klegend({'no drive','sweet spot drive'},'Location','southwest');

end

% Flux susceptibility of the dressed |g,0> to |g,1> cavity transition, Eq. (4.19)
function dnm=suscept(H_hilb,N_hilb,D_hilb,omega0,idx0,idx1)
dw=2*pi*1e3; shifts=[-dw dw]; frqs=zeros(1,2);
for k=1:2
    [V,E]=eig(full(H_hilb+omega0*D_hilb+shifts(k)*N_hilb));
    [~,k0]=max(abs(V(idx0,:))); [~,k1]=max(abs(V(idx1,:)));
    frqs(k)=E(k1,k1)-E(k0,k0);
end
dnm=(frqs(2)-frqs(1))/(2*dw);
end

