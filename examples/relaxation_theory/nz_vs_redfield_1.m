% Nakajima-Zwanzig relaxation theory against Redfield theory for a
% two-spin system with dipolar and CSA cross-correlations. The three
% superoperators compared are the off-shell NZ kernel (resolvent form),
% the on-shell NZ kernel (back-rotated form), and Redfield theory. The
% on-shell kernel at zero shift reproduces Redfield theory exactly; the
% off-shell kernel agrees with Redfield theory on the zero-frequency
% subspace of the coherent Liouvillian and differs in first order in
% omega*tau_c on coherences; a lifetime shift suppresses all rates by
% pushing the kernel off the real axis.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il
% talos@spindynamics.org

function nz_vs_redfield_1()

% Magnet and isotopes
sys.magnet=14.1;
sys.isotopes={'1H','13C'};

% Chemical shielding tensors
inter.zeeman.eigs={[7  15 -22]
                   [11 18 -29]};
inter.zeeman.euler={[pi/3 pi/4 pi/5]
                    [pi/6 pi/7 pi/8]};

% Coordinates for the dipolar coupling
inter.coordinates={[0.00 0.00 0.00]
                   [0.00 0.00 1.02]};

% Relaxation theory common settings
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.rlx_dfs='keep';

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Redfield superoperator
inter_rf=inter; inter_rf.relaxation={'redfield'}; inter_rf.tau_c={200e-12};
R_rf=full(relaxation(basis(create(sys,inter_rf),bas)));

% On-shell NZ superoperator at zero shift
inter_on=inter; inter_on.relaxation={'naka-zwan'}; inter_on.tau_c={200e-12};
inter_on.nz_onshell=true; inter_on.nz_shift=0;
R_on=full(relaxation(basis(create(sys,inter_on),bas)));

% Off-shell NZ superoperator at zero shift
inter_off=inter_on; inter_off.nz_onshell=false;
spin_system=basis(create(sys,inter_off),bas);
R_off=full(relaxation(spin_system));

% On-shell kernel at zero shift is Redfield theory
disp(['on-shell NZ vs Redfield, relative difference: ' ...
      num2str(norm(R_on-R_rf,1)/norm(R_rf,1))]);

% Off-shell kernel agrees with Redfield theory on the kite
V0=null(full(hamiltonian(assume(spin_system,'labframe'))));
disp(['off-shell NZ vs Redfield on the zero-frequency subspace: ' ...
      num2str(norm(R_off*V0-R_rf*V0,1)/norm(R_rf*V0,1))]);
disp(['off-shell NZ vs Redfield, whole superoperator: ' ...
      num2str(norm(R_off-R_rf,1)/norm(R_rf,1))]);

% Off-shell difference is first order in omega*tau_c
tau_grid=[2.5e-12 5e-12 10e-12 20e-12]; rel_diff=zeros(1,4);
for n=1:4
    inter_a=inter_rf; inter_a.tau_c={tau_grid(n)};
    R_a=full(relaxation(basis(create(sys,inter_a),bas)));
    inter_b=inter_off; inter_b.tau_c={tau_grid(n)};
    R_b=full(relaxation(basis(create(sys,inter_b),bas)));
    rel_diff(n)=norm(R_b-R_a,1)/norm(R_a,1);
end
disp('off-shell NZ vs Redfield difference against tau_c:');
disp([tau_grid; rel_diff]);

% Lifetime shift suppresses relaxation rates
shift_grid=[0 1e9 1e10 1e11]; max_rates=zeros(1,4);
for n=1:4
    inter_c=inter_off; inter_c.nz_shift=shift_grid(n);
    R_c=full(relaxation(basis(create(sys,inter_c),bas)));
    max_rates(n)=max(abs(diag(R_c)));
end
disp('largest relaxation rate against the lifetime shift, Hz:');
disp([shift_grid; max_rates]);

% Plot the two trends
kfigure(); scale_figure([1.75 0.75]);
subplot(1,2,1); plot(tau_grid,rel_diff,'-o');
kxlabel('$\tau_{c}$, seconds'); kgrid;
kylabel('rel. difference from Redfield');
subplot(1,2,2); semilogx(shift_grid(2:end),max_rates(2:end),'-o');
kxlabel('lifetime shift, Hz'); kgrid;
kylabel('largest relaxation rate, Hz');

end

