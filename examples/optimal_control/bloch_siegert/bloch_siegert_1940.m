% The Bloch-Siegert resonance displacement experiment of 1940 in
% silico. A single proton is driven by a linearly polarised field
% whose counter-rotating component displaces the resonance: at the
% transfer peak, the carrier frequency magnitude exceeds the Larmor
% frequency magnitude by w1^2/4w0. The displacement is measured in
% three independent ways: as the transfer peak of the corrected
% optimal control module, and as the nutation frequency minimum of
% the exact laboratory dynamics - on resonance the effective field
% in the rotating frame is smallest, and with it the nutation
% frequency, which is read off the rotation angle of the propagator
% accumulated over one carrier period - swept both in the spin
% frequency and in the carrier frequency. A final section measures
% the accuracy of the second-order effective Hamiltonian against
% the exact nutation frequencies: the residual follows the
% w1^3/16w0^2 law exactly at resonance and is fourth order in the
% field off resonance.
%
% A cautionary note: a transfer curve measured at weak drive does
% not locate the displacement correctly, because the ripple that
% the counter-rotating field puts on the magnetisation trajectory
% exceeds the secular variation of the transfer; the nutation
% frequency is the correct observable on the exact side, and the
% corrected module agrees with it.
%
% Calculation time: minutes.
%
% ilya.kuprov@weizmann.ac.il

function bloch_siegert_1940()

% Larmor frequency 2*pi*1e6 rad/s via the magnet value
sys.magnet=2*pi*1e6/spin('1H');

% Spin system
sys.isotopes={'1H'};

% Chemical shifts, ppm
inter.zeeman.scalar={0.0};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Drive amplitude and predicted displacement
w_s=spin_system.inter.basefrqs(1);
u=0.02*abs(w_s); d_shift=u^2/(4*abs(w_s));

% Offset grid around the predicted displacement
offs=d_shift*(-1:0.5:3); t_half=pi/u;

% Control operators and states
Lxo=operator(spin_system,'Lx','1H');
Lyo=operator(spin_system,'Ly','1H');
Lzo=operator(spin_system,'Lz','1H');
rho_z=state(spin_system,{'Lz'},{1});
rho_z=rho_z/norm(full(rho_z),2);

% Transfer sweep through the corrected module
q_mod=zeros(size(offs));
for m=1:numel(offs)

    % Control structure at the current spin offset
    ctrl=struct();
    ctrl.drifts={{offs(m)*Lzo}};
    ctrl.operators={Lxo,Lyo};
    ctrl.isotopes={'1H'}; ctrl.channels=[1;1];
    ctrl.rho_init={rho_z}; ctrl.rho_targ={rho_z};
    ctrl.pwr_levels=u; ctrl.pulse_dt=(t_half/100)*ones(1,100);
    ctrl.method='lbfgs'; ctrl.fidelity='real';
    ctrl.integrator='rectangle';
    ctrl.bsiegert=true();

    % Transfer after a half nutation period of constant drive
    ssm=optimcon(spin_system,ctrl);
    ssm.control.return_traj=true();
    [traj,~]=ensemble([ones(1,100); zeros(1,100)],ssm);
    q_mod(m)=1-real(rho_z'*traj{1}.forward(:,end));

end
pk_mod=parab_peak(offs,q_mod);
disp(['module transfer peak at ' num2str(pk_mod/d_shift,'%.4f') ...
      ' shift units (analytic: +1)']);

% Single-spin Zeeman and transverse operators in Hilbert space
bas_lab.formalism='zeeman-hilb'; bas_lab.approximation='none';
ss_lab=create(sys,inter); ss_lab=basis(ss_lab,bas_lab);
Sx=full(operator(ss_lab,'Lx','1H'));
Sz=full(operator(ss_lab,'Lz','1H'));

% Exact nutation frequency minimum, spin frequency swept
offs_sc=(0.02^2/4)*(-2:0.5:2); om_sc=zeros(size(offs_sc));
for m=1:numel(offs_sc)
    om_sc(m)=nut_freq(-1+offs_sc(m),-1,0.02,Sx,Sz);
end
[~,idx]=min(om_sc); idx=min(max(idx,2),numel(om_sc)-1);
cf=polyfit(offs_sc(idx-1:idx+1),om_sc(idx-1:idx+1),2);
pk_spin=-cf(2)/(2*cf(1));
disp(['exact nutation frequency minimum, spin-swept: ' ...
      num2str(pk_spin/(0.02^2/4),'%.4f') ' shift units (analytic: +1)']);

% Exact nutation frequency minimum, carrier frequency swept
om_c=zeros(size(offs_sc));
for m=1:numel(offs_sc)
    om_c(m)=nut_freq(-1,-1+offs_sc(m),0.02,Sx,Sz);
end
[~,idx]=min(om_c); idx=min(max(idx,2),numel(om_c)-1);
cf=polyfit(offs_sc(idx-1:idx+1),om_c(idx-1:idx+1),2);
pk_car=-cf(2)/(2*cf(1));
disp(['exact nutation frequency minimum, carrier-swept: ' ...
      num2str(pk_car/(0.02^2/4),'%.4f') ' shift units (analytic: -1)']);

% Residual of the effective Hamiltonian on resonance
rats=[0.04 0.02 0.01]; resid=zeros(1,3);
for m=1:3
    om_ex=nut_freq(-1,-1,rats(m),Sx,Sz);
    delta=rats(m)^2/(2*(-2));
    resid(m)=abs(om_ex-sqrt(delta^2+rats(m)^2));
    disp(['on resonance, r=' num2str(rats(m)) ...
          ': residual ' num2str(resid(m),'%.3e') ...
          ', cubic law ' num2str(rats(m)^3/16,'%.3e')]);
end

% Residual of the effective Hamiltonian off resonance
resid_o=zeros(1,2);
for m=1:2
    om_ex=nut_freq(-1,-0.5,rats(m+1),Sx,Sz);
    delta=(-1+0.5)+rats(m+1)^2/(2*(-1.5));
    om_mod=sqrt(delta^2+rats(m+1)^2);
    pred=mod(om_mod*(2*pi/0.5),2*pi); pred=min(pred,2*pi-pred);
    resid_o(m)=abs(om_ex-pred/(2*pi/0.5));
end
disp(['off resonance residual scaling per halving: ' ...
      num2str(resid_o(1)/resid_o(2),'%.1f') ' (fourth order: 16)']);

% Check the displacement measurements
if abs(pk_mod/d_shift-1)>0.1
    error('module displacement disagrees with analytics.');
end
if abs(pk_spin/(0.02^2/4)-1)>0.1
    error('exact spin-swept displacement disagrees with analytics.');
end
if abs(pk_car/(0.02^2/4)+1)>0.1
    error('exact carrier-swept displacement disagrees with analytics.');
end

% Check the residual laws
if max(abs(resid./(rats.^3/16)-1))>0.05
    error('on-resonance residual does not follow the cubic law.');
end
if (resid_o(1)/resid_o(2)<14)||(resid_o(1)/resid_o(2)>18)
    error('off-resonance residual is not fourth order.');
end
disp('Bloch-Siegert 1940 displacement confirmed in all three measurements');

end

% Nutation frequency from the propagator over one carrier period
function om=nut_freq(w_sp,w_cv,u,Sx,Sz)
t_car=2*pi/abs(w_cv); nstep=16384; dt=t_car/nstep;
U=eye(2);
for k=1:nstep
    t_mid=(k-0.5)*dt;
    H=w_sp*Sz+2*u*cos(w_cv*t_mid)*Sx;
    U=expm(-1i*H*dt)*U;
end
phs=angle(eig(U)); rot_ang=abs(phs(1)-phs(2));
rot_ang=min(rot_ang,2*pi-rot_ang);
om=rot_ang/t_car;
end

% Parabolic peak position from the three samples around the extremum
function pk=parab_peak(x,q)
[~,idx]=max(q); idx=min(max(idx,2),numel(q)-1);
cf=polyfit(x(idx-1:idx+1),q(idx-1:idx+1),2);
pk=-cf(2)/(2*cf(1));
end

