% Optimal control pulse optimisation with Bloch-Siegert shift corrections
% switched on. A single proton with a Larmor frequency of 1 MHz is driven
% at a significant fraction of its Larmor frequency -- a regime where the
% counter-rotating component of the control field shifts the resonance ap-
% preciably. A 90-degree pulse is optimised using LBFGS-GRAPE algorithm. 
% For comparison, the same pulse is also optimised with the corrections
% switched off and then evaluated in the corrected model.
%
% Note: the offset ensemble is deliberately applied through a transverse
% Lx operator rather than the usual Lz. The resulting static transverse
% fields of up to half the drive amplitude stress the Bloch-Siegert
% corrections much harder than resonance offsets would: the fidelity gap
% between the corrected and the uncorrected pulse is about six times
% larger than with an Lz offset ensemble.
%
% Calculation time: minutes.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=features_bss.m>

function features_bss()

% Larmor frequency of 1 MHz
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

% Deliberately transverse (Lx) offset ensemble, see the header note
control.off_ops={operator(spin_system,'Lx','1H')};
control.offsets={linspace(-1e5,1e5,11)};
control.isotopes={'1H'};
control.channels=[1;1];
control.rho_init={rho_init};
control.rho_targ={rho_targ};
control.pwr_levels=0.2*abs(spin_system.inter.basefrqs(1));
control.pulse_dt=(8.0*pi/control.pwr_levels/50)*ones(1,50);
control.method='lbfgs';
control.fidelity='real';
control.max_iter=100;

% Deterministic initial guess
guess=[linspace(0.1,0.5,50); 0.05*ones(1,50)];

% Optimisation
control.bsiegert=true();
spin_system=optimcon(spin_system,control);
pulse=fmaxnewton(spin_system,@grape_xy,guess);

% Corrected pulse in the corrected model
[~,fid_corr]=ensemble(pulse,spin_system);

% Optimisation
control.bsiegert=false();
spin_system=optimcon(spin_system,control);
pulse_off=fmaxnewton(spin_system,@grape_xy,guess);

% Uncorrected pulse in the corrected model
control.bsiegert=true();
spin_system=optimcon(spin_system,control);
[~,fid_off]=ensemble(pulse_off,spin_system);

% Report the results
report(spin_system,'---------------');
report(spin_system,['corrected pulse, corrected model:   ' ...
                     num2str(fid_corr,'%.6f')]);
report(spin_system,['uncorrected pulse, corrected model: ' ...
                     num2str(fid_off,'%.6f')]);

end

