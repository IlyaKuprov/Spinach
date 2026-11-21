% Switch a heteronuclear scalar interaction off and back on.
% The proton coherence accumulates conditional phase only while J is on.
%
% talos@spindynamics.org

function interaction_switch()

% Specify a proton and a carbon with an 80 Hz scalar coupling
sys.magnet=9.4;
sys.isotopes={'1H','13C'};
sys.parallel={'processes',1};
sys.output='hush';
inter.zeeman.scalar={0,0};
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=80;
bas.formalism='zeeman-hilb';
bas.approximation={'none'};

% Precompute the two rotating-frame generators before any evolution
spin_system=basis(create(sys,inter),bas);
H=hamiltonian(assume(spin_system,'nmr'));
inter.coupling.scalar{1,2}=0;
uncoupled=basis(create(sys,inter),bas);
H_off=hamiltonian(assume(uncoupled,'nmr'));

% Prepare proton coherence and propagate through the three intervals
rho=full(state(spin_system,'L+','1H'));
coil=full(coil_state(spin_system,'L+','1H'));
dt=1e-4;
nsteps=[20 30 20];
generators={H,H_off,H};
signal=complex(zeros(1,1+sum(nsteps)));
signal(1)=trace(coil'*rho);
first=1;
for n=1:numel(nsteps)

    % Reuse one propagator for all samples in this interval
    P=expm(-1i*full(generators{n})*dt);
    for k=1:nsteps(n)
        rho=P*rho*P';
        signal(first+k)=trace(coil'*rho);
    end
    first=first+nsteps(n);

end

% Compare all samples to the conditional phase accumulated during J-on time
on_time=dt*[0 cumsum([ones(1,20) zeros(1,30) ones(1,20)])];
expected=cos(pi*80*on_time);
signal=signal/signal(1);
assert(max(abs(signal-expected))<1e-10);

% Plot both the held coherence and its resumed scalar evolution
kfigure();
plot((0:sum(nsteps))*dt,real(signal),(0:sum(nsteps))*dt,expected,'--');
kxlabel('time, seconds');
kylabel('normalised proton coherence');
klegend({'Switched generators','Analytic evolution'});

end


