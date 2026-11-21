% One-spin state transfer with circular peak RF and a pulse energy cap.
% A frozen sample is included in both feasibility checks and energy.
%
% talos@spindynamics.org

function peak_rms_design()

% Build a full Hilbert-space control problem with two RF quadratures
sys.magnet=0;
sys.isotopes={'1H'};
sys.parallel={'processes',1};
sys.output='hush';
inter.zeeman.scalar={0};
bas.formalism='zeeman-hilb';
bas.approximation={'none'};
spin_system=basis(create(sys,inter),bas);
S=pauli(2);
control.isotopes={'1H'};
control.channels=[1 1];
control.operators={S.x,S.y};
control.rho_init={S.z};
control.rho_targ={S.x};
control.drifts={{sparse(2,2)}};
control.pwr_levels=1;
control.pulse_dt=[0.05 0.15];
control.method='lbfgs';
control.penalties={'none'};
control.p_weights=0;
control.plotting={};
spin_system=optimcon(spin_system,control);

% Solve the bounded design and check the known resonant optimum
peak=2;
fixed=[0;0];
energy=control.pulse_dt(2)*peak^2;
[free,fidelity]=fit_pulse(spin_system,peak,energy,fixed);
assert(norm(free-[0;peak])<1e-4);
fprintf('Peak %.6f rad/s, energy %.6f rad^2/s, fidelity %.6f\n',...
        norm(free),control.pulse_dt(2)*sum(free.^2),fidelity);

% Check an energy-limited solution starting from zero RF
[free,~]=fit_pulse(spin_system,peak,energy/4,fixed);
assert(norm(free-[0;peak/2])<1e-4);

% Refuse a frozen waveform that already exceeds the energy budget
try
    fit_pulse(spin_system,peak,energy/8,[0;peak]);
    refused=false;
catch err
    refused=strcmp(err.identifier,'Spinach:peak_rms:infeasible');
end
assert(refused);

end

% Solve the constrained two-quadrature design without clipping trial points
function [free,fidelity]=fit_pulse(spin_system,peak,energy,fixed)

% Deduct the frozen sample energy before optimising the free quadratures
dt=spin_system.control.pulse_dt;
remaining=energy-dt(1)*sum(fixed.^2);
if norm(fixed)>peak||remaining<0
    error('Spinach:peak_rms:infeasible','frozen RF violates the peak or energy limit.');
end

% Use the explicit nonlinear constraints with a feasible zero start
options=optimoptions('fmincon','Display','off','Algorithm','sqp');
objective=@(free)physical_cost(free,fixed,spin_system);
constraint=@(free)limits(free,peak,remaining,dt(2));
[free,cost,exitflag]=fmincon(objective,zeros(2,1),[],[],[],[],[],[],constraint,options);
assert(exitflag>0&&all(limits(free,peak,remaining,dt(2))<=options.ConstraintTolerance));
fidelity=-cost;

end

% Minimise the negative physical GRAPE fidelity
function cost=physical_cost(free,fixed,spin_system)
[~,fidelity]=grape_xy([fixed free],spin_system);
cost=-fidelity(1);
end

% Enforce the circular peak and energy constraints
function [ineq,eq]=limits(free,peak,remaining,dt)
norm_sq=sum(free.^2);
ineq=[norm_sq-peak^2; dt*norm_sq-remaining];
eq=[];
end


