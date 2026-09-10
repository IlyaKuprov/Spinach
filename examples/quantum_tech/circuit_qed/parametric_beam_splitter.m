% Parametrically activated beam-splitter interaction between two
% cavities coupled through a flux-modulated transmon, Sec. 2.3.1
% of Yunwei Lu, PhD thesis, Northwestern University (2026). The
% linear part of the Hamiltonian is diagonalised into dressed
% normal modes; in that basis the transmon frequency modulation
% delta_omega_b*cos(omega_d*t+phi_d)*b'b acquires an off-diagonal
% term that swaps a photon between the two cavity-like modes when
% omega_d matches their dressed frequency difference. The swap
% rate G_ac=delta_omega_b*abs(eta_ac)/2 of Eq. (2.69) is checked
% against the population dynamics of a single photon placed into
% the dressed cavity mode a. The simulation runs in the frame
% rotating at the bare transmon frequency; the propagator over one
% drive period is assembled from piecewise-constant steps and then
% applied stroboscopically. Without the modulation, the photon
% stays in the dressed mode it was placed into.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function parametric_beam_splitter()

% Magnet field
sys.magnet=0;

% Cavity a, transmon b, and cavity c
sys.isotopes={'C3','T3','C3'};

% Mode frequencies as detunings from the transmon frequency
inter.modes.carriers={5.5e9 5.5e9 5.5e9};
inter.modes.frqs={0.5e9 0 1.5e9};
inter.modes.anharms={[] -200e6 []};

% Static exchange couplings of the transmon to the cavities
inter.modes.exchange=cell(3,3);
inter.modes.exchange{1,2}=50e6;
inter.modes.exchange{2,3}=50e6;

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Static Hamiltonian and transmon number operator, Eqs. (2.62) and (2.63)
H0=hamiltonian(assume(spin_system,'cavity'));
Nb=operator(spin_system,'N',2);

% Indices of the bare single-excitation states
idx=zeros(1,3);
for k=1:3
    levels={'BL1','BL1','BL1'}; levels{k}='BL2';
    idx(k)=find(diag(state(spin_system,levels,{1,2,3}))>0.5);
end

% Dressed normal modes from the single-excitation block, Eq. (2.65)
[vecs,vals]=eig(full(H0(idx,idx))); frq_dressed=diag(vals);

% Label the dressed modes by their dominant bare mode
[~,dominant]=max(abs(vecs),[],1); [~,order]=sort(dominant);
vecs=vecs(:,order); frq_dressed=frq_dressed(order);

% Transmon weights in the dressed modes, Eq. (2.66)
u_b=vecs(2,:);

% Modulation amplitude, phase, and drive frequency, Eq. (2.64)
mod_amp=2*pi*100e6; phi_d=pi/3;
frq_drive=frq_dressed(3)-frq_dressed(1);

% Beam-splitter rate from Eq. (2.69)
eta_ac=conj(u_b(1))*u_b(3); rate_theory=mod_amp*abs(eta_ac)/2;

% Free evolution propagator over one drive period
period=2*pi/frq_drive; P0=propagator(spin_system,H0,period);

% Piecewise-constant driven propagator over one drive period
nsteps=50; dt=period/nsteps; P1=speye(size(H0,1));
for k=1:nsteps
    t_mid=(k-0.5)*dt;
    H=H0+mod_amp*cos(frq_drive*t_mid+phi_d)*Nb;
    P1=propagator(spin_system,H,dt)*P1;
end

% Stroboscopic trajectories covering a full swap, free and driven
nper=ceil(1.1*pi/(rate_theory*period)); pops=zeros(3,nper+1,2);
props={P0 P1};
for m=1:2
    psi=zeros(size(H0,1),1); psi(idx)=vecs(:,1);
    pops(:,1,m)=abs(vecs'*psi(idx)).^2;
    for n=1:nper
        psi=props{m}*psi; pops(:,n+1,m)=abs(vecs'*psi(idx)).^2;
    end
end

% Time axis and driven population of dressed mode c
time_axis=period*(0:nper); pop_c=pops(3,:,2);

% Validate the completeness of the driven photon transfer
if max(pop_c)<0.9
    error('driven photon transfer between the cavities is incomplete.');
end

% Interval during which the mode c population exceeds one half
n_up=find(pop_c>0.5,1); n_dn=n_up+find(pop_c(n_up:end)<0.5,1)-1;
if isempty(n_dn)
    error('mode c population does not return below one half within the simulated window.');
end
t_up=interp1(pop_c(n_up-1:n_up),time_axis(n_up-1:n_up),0.5);
t_dn=interp1(pop_c(n_dn-1:n_dn),time_axis(n_dn-1:n_dn),0.5);

% Population sin(G*t)^2 stays above one half for pi/(2G)
rate_numer=pi/(2*(t_dn-t_up));
disp(['swap rate, numerical:  ' num2str(rate_numer/(2*pi*1e3)) ' kHz']);
disp(['swap rate, Eq. (2.69): ' num2str(rate_theory/(2*pi*1e3)) ' kHz']);

% Validate the swap rate against Eq. (2.69)
if abs(rate_numer/rate_theory-1)>0.1
    error('numerical swap rate deviates from Eq. (2.69).');
end

% Validate that the photon stays put without the drive
if min(pops(1,:,1))<0.999
    error('photon leaves the dressed cavity mode without the drive.');
end

% Validate confinement to the single-excitation manifold
if max(abs(sum(pops,1)-1),[],'all')>1e-6
    error('dressed single-excitation populations are not conserved.');
end

% Plot dressed mode populations with and without the drive
kfigure(); scale_figure([2.0 0.75]);
subplot(1,2,1); plot(1e6*time_axis,pops(:,:,1)','LineWidth',1.5);
axis tight; ylim([0 1]); kgrid; kxlabel('time, $\mu$s'); kylabel('dressed mode population');
ktitle('no modulation');
klegend({'mode a','mode b','mode c'},'Location','Best');
subplot(1,2,2); plot(1e6*time_axis,pops(:,:,2)','LineWidth',1.5);
axis tight; ylim([0 1]); kgrid; kxlabel('time, $\mu$s'); kylabel('dressed mode population');
ktitle('parametric modulation at $\tilde{\omega}_c-\tilde{\omega}_a$');
klegend({'mode a','mode b','mode c'},'Location','Best');

end

