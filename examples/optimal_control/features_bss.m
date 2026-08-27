% Optimal control pulse optimisation with Bloch-Siegert shift corrections
% switched on. A single proton with a Larmor frequency of 1 MHz is driven
% at a tenth of its Larmor frequency - a regime where the counter-rotating
% component of the control field shifts the resonance appreciably. A 90-
% degree pulse is optimised with the LBFGS-GRAPE algorithm in a model that
% carries the Bloch-Siegert corrections through the response operators
% that optimcon.m builds using bss_ops.m from the control channel map
% when control.bsiegert is set; the corrections enter each time
% slice through the squares of the control amplitudes, and the fidelity
% gradient is pulled back through that quadratic map inside the GRAPE
% engines.
%
% For comparison, the same pulse is also optimised with the corrections
% switched off and then evaluated in the corrected model: the fidelity
% penalty of ignoring the Bloch-Siegert physics is printed.
%
% Calculation time: minutes.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=features_bss.m>

function features_bss()

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

% Set up and normalise the initial state
rho_init=state(spin_system,{'Lz'},{1});
rho_init=rho_init/norm(full(rho_init),2);

% Set up and normalise the target state
rho_targ=state(spin_system,{'Lx'},{1});
rho_targ=rho_targ/norm(full(rho_targ),2);

% Drift Liouvillian is zero on resonance
dim=size(spin_system.bas.basis,1);

% Control parameters
control.drifts={{sparse(dim,dim)}};
control.operators={operator(spin_system,'Lx','1H'),...
                   operator(spin_system,'Ly','1H')};
control.isotopes={'1H'};
control.channels=[1;1];
control.rho_init={rho_init};
control.rho_targ={rho_targ};
control.pwr_levels=0.1*abs(spin_system.inter.basefrqs(1));
control.pulse_dt=(2.5*pi/control.pwr_levels/50)*ones(1,50);
control.method='lbfgs';
control.fidelity='real';
control.max_iter=100;

% Switch the Bloch-Siegert corrections on
control.bsiegert=true();

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Deterministic initial guess
guess=[linspace(0.1,0.5,50); 0.05*ones(1,50)];

% Optimise with Bloch-Siegert corrections on
pulse=fmaxnewton(spin_system,@grape_xy,guess);

% Fidelity of the corrected pulse in the corrected model
[~,fid_corr]=ensemble(pulse,spin_system);
report(spin_system,['corrected pulse, corrected model:   ' ...
                    num2str(fid_corr,'%.6f')]);

% Rebuild the problem with the corrections off
control.bsiegert=false();
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system=optimcon(spin_system,control);

% Optimise with Bloch-Siegert corrections off
pulse_off=fmaxnewton(spin_system,@grape_xy,guess);

% Evaluate the uncorrected pulse in the corrected model
control.bsiegert=true();
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system=optimcon(spin_system,control);
[~,fid_off]=ensemble(pulse_off,spin_system);
report(spin_system,['uncorrected pulse, corrected model: ' ...
                    num2str(fid_off,'%.6f')]);

end

