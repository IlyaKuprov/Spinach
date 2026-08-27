% Ramsey shifts of other spins under an off-resonant drive. A proton
% channel drives a three-spin system in which 13C and 15N are far off
% resonance: each of their Zeeman frequencies shifts by
%
%             delta_n=w_n*w1n^2/(w_n^2-w_c^2)
%
% where w_n is the signed Zeeman frequency of the nucleus, w_c is the
% signed carrier frequency, and w1n is the drive amplitude scaled by
% the ratio of the magnetogyric ratios. The negative magnetogyric
% ratio of 15N makes its shift precess the opposite way to 13C. The
% accumulated phases of the off-resonant nuclei are compared with the
% analytic formula, and the quadratic scaling of the shift in the
% drive amplitude and its inverse scaling in the magnet field are
% demonstrated.
%
% Calculation time: seconds.
%
% ilya.kuprov@weizmann.ac.il

function ramsey_shifts()

% Off-resonant nucleus phases at the base amplitude and field
[ph_c,ph_n,ph_ca,ph_na]=shift_phases(14.1,2*pi*25e3);
disp(['13C phase: ' num2str(ph_c,'%.9f') ...
      ' rad, analytic ' num2str(ph_ca,'%.9f') ' rad']);
disp(['15N phase: ' num2str(ph_n,'%.9f') ...
      ' rad, analytic ' num2str(ph_na,'%.9f') ' rad']);

% Check the phases and the sign flip
if abs(ph_c-ph_ca)/abs(ph_ca)>1e-6
    error('13C phase disagrees with the Ramsey formula.');
end
if abs(ph_n-ph_na)/abs(ph_na)>1e-6
    error('15N phase disagrees with the Ramsey formula.');
end
if sign(ph_c)==sign(ph_n)
    error('15N shift must precess opposite to 13C.');
end
disp('phases match the Ramsey formula, sign flip confirmed');

% Quadratic scaling in the drive amplitude
ph_c2=shift_phases(14.1,2*pi*50e3);
disp(['amplitude doubling phase ratio: ' ...
      num2str(ph_c2/ph_c,'%.8f') ' (analytic: 4)']);
if abs(ph_c2/ph_c-4)>1e-5
    error('Ramsey shift is not quadratic in the drive amplitude.');
end

% Inverse scaling in the magnet field
ph_ch=shift_phases(14.1/2,2*pi*25e3);
disp(['field halving phase ratio: ' ...
      num2str(ph_ch/ph_c,'%.8f') ' (analytic: 2)']);
if abs(ph_ch/ph_c-2)>1e-5
    error('Ramsey shift is not inverse in the magnet field.');
end
disp('Ramsey shift physics confirmed');

end

% Phases accumulated by the off-resonant nuclei under a proton drive
function [ph_c,ph_n,ph_ca,ph_na]=shift_phases(b_field,amp)

% Three-spin system with two off-resonant isotopes
sys.magnet=b_field; sys.isotopes={'1H','13C','15N'};
sys.output='hush';
inter.zeeman.scalar={0.0,0.0,0.0};
bas.formalism='sphten-liouv'; bas.approximation='none';
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
dim=size(spin_system.bas.basis,1);

% Transverse magnetisation of the off-resonant nuclei
rho_c=state(spin_system,{'L+'},{2});
rho_c=rho_c/norm(full(rho_c),2);
rho_n=state(spin_system,{'L+'},{3});
rho_n=rho_n/norm(full(rho_n),2);
rho_init=(rho_c+rho_n)/norm(full(rho_c+rho_n),2);

% Constant proton drive for twenty milliseconds
t_tot=0.02;
ctrl.drifts={{sparse(dim,dim)}};
ctrl.operators={operator(spin_system,'Lx','1H'),...
                operator(spin_system,'Ly','1H')};
ctrl.isotopes={'1H'}; ctrl.channels=[1;1];
ctrl.rho_init={rho_init}; ctrl.rho_targ={rho_init};
ctrl.pwr_levels=amp; ctrl.pulse_dt=(t_tot/50)*ones(1,50);
ctrl.method='lbfgs'; ctrl.fidelity='real';
ctrl.integrator='rectangle';
ctrl.bsiegert=true();
spin_system=optimcon(spin_system,ctrl);
spin_system.control.return_traj=true();
[traj,~]=ensemble([ones(1,50); zeros(1,50)],spin_system);

% Accumulated phases of the off-resonant nuclei
ph_c=angle(rho_c'*traj{1}.forward(:,end));
ph_n=angle(rho_n'*traj{1}.forward(:,end));

% Analytic Ramsey phases
gam=spin_system.inter.gammas;
bfrq=spin_system.inter.basefrqs; wc=bfrq(1);
delta_c=(gam(2)/gam(1))^2*amp^2*bfrq(2)/(bfrq(2)^2-wc^2);
delta_n=(gam(3)/gam(1))^2*amp^2*bfrq(3)/(bfrq(3)^2-wc^2);
ph_ca=angle(exp(-1i*delta_c*t_tot));
ph_na=angle(exp(-1i*delta_n*t_tot));

end

