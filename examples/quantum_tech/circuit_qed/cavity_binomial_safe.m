% Binomial bosonic code |0L>=(|0>+|4>)/sqrt(2), |1L>=|2> in a cavity
% dispersively coupled to a flux-tunable transmon ancilla, and the
% protection of its coherences from 1/f flux noise by a Stark-assis-
% ted flux-noise evasion (SAFE) drive on the transmon, Sec. 4.4.1 and
% Fig. 4.4(a,b) of Yunwei Lu's PhD thesis (Northwestern University,
% 2026). The flux noise dephasing rates of the code and error space
% coherences are computed from the flux sensitivities of the dressed
% cavity transition frequencies as functions of the transmon-drive
% detuning, Eq. (4.93); at the common minimum the logical state |+L>
% is then propagated for 300 microseconds along 1/f flux noise tra-
% jectories under the Lindblad master equation, with and without the
% drive, and the decoherence-only infidelity of Eq. (4.94) and the
% Wigner function of the cavity state are reported.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function cavity_binomial_safe()

% Transmon anharmonicity placing the sweet spot of Eq. (C.22) at -30 MHz for a 10 MHz drive, Hz
anharm=-67e6;

% Transmon-cavity exchange coupling and detuning giving a 0.5 MHz dispersive shift, Hz
g_bc=86e6; delta_bc=1.414e9;

% Dispersive shift and the sensitivities of the cavity terms to the transmon frequency
chi=2*(g_bc/delta_bc)^2*anharm; sens_c=(g_bc/delta_bc)^2; sens_x=-4*anharm*g_bc^2/delta_bc^3;

% Transmon frequency sensitivity to flux, rad/s per flux quantum
dwb_dphi=2*pi*6e9;

% Flux noise amplitude in flux quanta and the ultraviolet cutoff, Hz
noise_amp=1e-5; f_uv=5e7;

% Drive amplitude, rad/s, and the transmon-drive detunings to scan, Hz
omega0=2*pi*10e6; detunings=-(20:1:80)*1e6;

% Transmon and cavity relaxation times, seconds
t1_b=50e-6; t1_c=20e-3;

% Time step, number of steps, infidelity sampling stride, and trajectory count
dt=1e-8; nsteps=30000; stride=300; ntraj=100;

% Length of the noise synthesis grid and the infrared cutoff at its frequency resolution, Hz
nlong=2^19; f_ir=1/(nlong*dt);

% Fock state pairs whose coherences matter for the code, Sec. 4.4.1
pairs=[2 0; 4 2; 4 0; 3 0; 4 3; 2 1; 3 1];

% Three-level transmon and a five-level cavity, both in their own frames
sys.magnet=0; sys.isotopes={'T3','C5'};
inter.modes.frqs={0 0};
inter.modes.anharms={anharm []};
inter.modes.lifetimes={t1_b t1_c};
inter.modes.kerr=cell(2,2); inter.modes.kerr{1,2}=chi;
inter.temperature=0;
bas.approximation='none';

% Hilbert space twin for the dressed states
bas.formalism='zeeman-hilb';
ss_hilb=create(sys,inter);
ss_hilb=basis(ss_hilb,bas);

% Hamiltonian without the detuning term, detuning, drive, and flux noise operators
H_hilb=hamiltonian(assume(ss_hilb,'labframe'));
num_b=operator(ss_hilb,'N',1);
drive_op=operator(ss_hilb,'C',1)+operator(ss_hilb,'A',1);
noise_op=num_b+sens_c*operator(ss_hilb,'N',2)+sens_x*operator(ss_hilb,{'N','N'},{1,2});

% Positions of the bare |g,n> states
idx=zeros(1,5);
for n=1:5
    idx(n)=find(diag(state(ss_hilb,{'BL1',['BL' int2str(n)]},{1,2}))>0.5);
end

% Flux noise dephasing rates of the code coherences with and without the drive, Eq. (4.93) with |ln(omega_ir*t)|=4
rates=zeros(size(pairs,1),numel(detunings)); dw=2*pi*1e3;
for n=1:numel(detunings)
    dressed=dressed_ens(H_hilb+2*pi*detunings(n)*num_b,noise_op,drive_op,omega0,idx,dw);
    rates(:,n)=noise_amp*dwb_dphi*sqrt(2*4)*abs(dressed(pairs(:,1)+1)-dressed(pairs(:,2)+1))';
end

% Operating detuning at the median of the rate minima, and the undriven rates there
[~,min_idx]=min(rates,[],2); op_idx=round(median(min_idx)); delta_bd=detunings(op_idx);
dressed=dressed_ens(H_hilb+2*pi*delta_bd*num_b,noise_op,drive_op,0,idx,dw);
rates_off=noise_amp*dwb_dphi*sqrt(2*4)*abs(dressed(pairs(:,1)+1)-dressed(pairs(:,2)+1))';

% Report the suppression at the operating detuning
disp(['rate minima at ' num2str(-detunings(min_idx)/1e6,'%.0f ') 'MHz, operating detuning ' num2str(-delta_bd/1e6) ' MHz']);
disp(['suppression factors at the operating detuning: ' num2str((rates_off./rates(:,op_idx))','%.1f ')]);

% Validate the sweet spot window
if any(min_idx==1)||any(min_idx==numel(detunings))
    error('dephasing rate minima are outside the detuning window.');
end
if any(rates_off./rates(:,op_idx)<5)
    error('SAFE drive does not suppress all code coherence dephasing rates fivefold.');
end

% Plot the dephasing rates, Fig. 4.4(a)
kfigure(); subplot(1,3,1); semilogy(-detunings/1e6,1e-9*rates'); kgrid;
kxlabel('$-\Delta_{bd}/2\pi$, MHz'); kylabel('dephasing rate, 1/ns'); ktitle('binomial code coherences');
klegend(cellstr(num2str(pairs,'$\\gamma_{%d,%d}$')),'Location','southeast');

% Liouville space for the dissipative dynamics
bas.formalism='zeeman-liouv';
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Drift, drive, and flux noise generators with the mode dissipators
L_drift=hamiltonian(assume(spin_system,'labframe'))+2*pi*delta_bd*operator(spin_system,'N',1)+...
        1i*relaxation(spin_system);
L_drive=operator(spin_system,'C',1)+operator(spin_system,'A',1);
L_noise=operator(spin_system,'N',1)+sens_c*operator(spin_system,'N',2)+...
        sens_x*operator(spin_system,{'N','N'},{1,2});

% Flux noise trajectories windowed from the long grid, transmon frequency units, rad/s
rng(1); dwb=zeros(ntraj,nsteps);
for m=1:ntraj
    dphi=pink_noise(noise_amp,dt,nlong,1,f_ir,f_uv); dwb(m,:)=dwb_dphi*dphi(1:nsteps);
end

% Frequency grid for the propagator table
ngrid=201; dw_max=max(abs(dwb(:))); dw_grid=linspace(-dw_max,dw_max,ngrid);
grid_idx=round((dwb+dw_max)*(ngrid-1)/(2*dw_max))+1;

% Time axis of the recorded infidelity
time_axis=dt*stride*(1:(nsteps/stride));

% Preallocate infidelities and final cavity states, and set the drive amplitudes of the two cases
infids=zeros(2,numel(time_axis)); rho_cav=zeros(5,5,3); drives=[0 omega0];

% Initial cavity state |+L> of the bare binomial code
psi_code=[1 0 sqrt(2) 0 1]'/2; rho_cav(:,:,1)=psi_code*psi_code';

% Loop over the undriven and the driven case
for k=1:2

    % Dressed states adiabatically connected to |g,0>, |g,2>, and |g,4>
    H_case=full(H_hilb+2*pi*delta_bd*num_b+drives(k)*drive_op); [vecs,vals]=eig(H_case);
    [~,c0]=max(abs(vecs(idx(1),:))); [~,c2]=max(abs(vecs(idx(3),:))); [~,c4]=max(abs(vecs(idx(5),:)));

    % Phases of the dressed states fixed by their bare state components
    cols=[c0 c2 c4]; bare=[idx(1) idx(3) idx(5)];
    for j=1:3
        vecs(:,cols(j))=vecs(:,cols(j))*abs(vecs(bare(j),cols(j)))/vecs(bare(j),cols(j));
    end

    % Logical |+L> state and its closed-system evolution
    psi=(vecs(:,c0)+vecs(:,c4))/2+vecs(:,c2)/sqrt(2); rho=hilb2liouv(psi*psi','statevec');
    psi_ideal=vecs*(exp(-1i*diag(vals)*time_axis).*(vecs'*psi));

    % Table of propagators over the frequency grid
    P=cell(1,ngrid);
    for n=1:ngrid
        P{n}=full(propagator(spin_system,L_drift+drives(k)*L_drive+dw_grid(n)*L_noise,dt));
    end

    % Ensemble average of the density matrix along the noise trajectories
    rho_avg=zeros(numel(rho),numel(time_axis));
    for m=1:ntraj
        rho_traj=rho;
        for n=1:nsteps
            rho_traj=P{grid_idx(m,n)}*rho_traj;
            if mod(n,stride)==0
                rho_avg(:,n/stride)=rho_avg(:,n/stride)+rho_traj;
            end
        end
    end
    rho_avg=rho_avg/ntraj;

    % Decoherence-only infidelity, Eq. (4.94)
    for n=1:numel(time_axis)
        rho_mat=reshape(rho_avg(:,n),size(H_case));
        infids(k,n)=1-sqrt(real(psi_ideal(:,n)'*rho_mat*psi_ideal(:,n)));
    end

    % Final cavity density matrix with the transmon traced out and the propagator truncation drift of the trace removed
    rho_end=reshape(rho_mat,[5 3 5 3]);
    for t=1:3
        rho_cav(:,:,k+1)=rho_cav(:,:,k+1)+squeeze(rho_end(:,t,:,t));
    end
    rho_cav(:,:,k+1)=rho_cav(:,:,k+1)/trace(rho_cav(:,:,k+1));

end

% Report and validate the infidelities
disp(['infidelity after ' num2str(1e6*time_axis(end)) ' us: without the drive ' num2str(infids(1,end),'%.3f') ...
      ', with the drive ' num2str(infids(2,end),'%.3f')]);
if infids(2,end)>infids(1,end)/3
    error('SAFE drive does not reduce the decoherence-induced infidelity threefold.');
end

% Wigner functions of the cavity states, Fig. 4.4(b), on a Fock basis padded to fit the displaced states
[re_grid,im_grid]=meshgrid(linspace(-3.5,3.5,71)); alpha=re_grid+1i*im_grid; wigners=zeros([size(alpha) 3]);
for k=1:3
    wigners(:,:,k)=wigner_fock(blkdiag(rho_cav(:,:,k),zeros(55)),alpha);
end

% Validate the Wigner function normalisation against the unit integral
if abs(sum(wigners(:,:,1),'all')*(re_grid(1,2)-re_grid(1,1))^2-1)>0.05
    error('Wigner function of the initial state does not integrate to unity.');
end

% Plot the infidelities and the Wigner functions
subplot(1,3,2); plot(1e6*time_axis,infids'); kgrid; kxlabel('time, $\mu$s'); kylabel('infidelity');
klegend({'no drive','sweet spot drive'},'Location','northwest'); ktitle('$|+_L\rangle$ decoherence');
subplot(1,3,3); contour(re_grid,im_grid,wigners(:,:,2),'LineWidth',1); hold on;
contour(re_grid,im_grid,wigners(:,:,3),'--','LineWidth',1); kgrid; axis square;
kxlabel('$I$'); kylabel('$Q$'); ktitle('Wigner function after 300 $\mu$s');
klegend({'no drive','sweet spot drive'},'Location','northeast');

end

% Flux susceptibilities of the dressed |g,n> energies, Eq. (4.19)
function suscept=dressed_ens(H_hilb,noise_op,drive_op,omega0,idx,dw)
energies=zeros(2,numel(idx)); shifts=[-dw dw];
for m=1:2
    [vecs,vals]=eig(full(H_hilb+omega0*drive_op+shifts(m)*noise_op));
    for j=1:numel(idx)
        [~,col]=max(abs(vecs(idx(j),:))); energies(m,j)=vals(col,col);
    end
end
suscept=(energies(2,:)-energies(1,:))/(2*dw);
end

