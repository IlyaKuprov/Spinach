% GRAPE preparation of a cavity Fock state through a dispersively
% coupled qubit using smooth band-limited drives. The controls are
% expanded in an orthonormal basis of slow sine and cosine waves,
% and the optimisation runs over the expansion coefficients, so
% that the resulting pulses are hardware-friendly smooth envelopes
% rather than free piecewise-constant switches. Compare with the
% piecewise-constant treatment in cavity_fock_grape_a.m; model and
% parameters follow the smooth pulse bosonic GRAPE example of the
% paraqeet package.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function cavity_fock_grape_b()

% Magnet field
sys.magnet=0;

% Truncated cavity mode and a qubit
sys.isotopes={'C4','E'};

% Cavity on resonance with its drive
inter.modes.frqs={0 []};

% Dispersive coupling to the qubit
inter.modes.dispersive=cell(2,2);
inter.modes.dispersive{1,2}=656.2e3;

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Drift Hamiltonian from the declared interactions
H=hamiltonian(assume(spin_system,'cavity'));

% Quadrature control operators of the cavity and the qubit
Cr=operator(spin_system,'C',1); An=operator(spin_system,'A',1);
Lx=operator(spin_system,'Lx',2); Ly=operator(spin_system,'Ly',2);
ops={(Cr+An)/2,1i*(Cr-An)/2,Lx,Ly};

% Vacuum cavity and cavity Fock state two, qubit in the upper level
rho_init=state(spin_system,{'BL1','ZL2'},{1,2});
rho_targ=state(spin_system,{'BL3','ZL2'},{1,2});
rho_init=rho_init/norm(rho_init,'fro');
rho_targ=rho_targ/norm(rho_targ,'fro');

% Band-limited orthonormal waveform basis
wave_set=[wave_basis('sine_waves',2,40) wave_basis('cosine_waves',3,40)]';

% Define control parameters
control.isotopes={'C4','E'};
control.channels=[1;1;2;2];
control.drifts={{H}};
control.operators=ops;
control.rho_init={rho_init};
control.rho_targ={rho_targ};
control.basis=wave_set;
control.pwr_levels=1.76828e7;
control.pulse_dt=33e-9*ones(1,40);
control.penalties={'NS'};
control.p_weights=0.001;
control.method='lbfgs';
control.max_iter=300;

% Plots during optimisation
control.plotting={'xy_controls'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Deterministic asymmetric coefficient guess
guess=[0.3 0.1 0.5 0.1 0.0; 0.0 0.1 0.0 0.1 0.0;
       0.3 0.0 0.5 0.0 0.1; 0.1 0.0 0.0 0.1 0.0];

% Run the optimisation, get basis coefficients
coeffs=fmaxnewton(spin_system,@grape_xy,guess);

% Reassemble the smooth time-domain pulse
pulse=coeffs*wave_set;

% Recompute the transfer fidelity by direct propagation
rho=rho_init;
for n=1:numel(spin_system.control.pulse_dt)
    slice_ham=H+spin_system.control.pwr_levels*...
                (pulse(1,n)*ops{1}+pulse(2,n)*ops{2}+...
                 pulse(3,n)*ops{3}+pulse(4,n)*ops{4});
    P=propagator(spin_system,slice_ham,spin_system.control.pulse_dt(n));
    rho=P*rho*P';
end
fid=real(trace(rho_targ'*rho));

% Validate the optimisation
if fid<0.80
    error('band-limited GRAPE optimisation did not converge.');
end

% Slice midpoint time axis of the optimised pulse
dts=spin_system.control.pulse_dt; time_axis=1e6*(cumsum(dts)-dts/2);

% Plot the smooth control envelopes
kfigure(); plot(time_axis,spin_system.control.pwr_levels*pulse'/(2*pi*1e6),'LineWidth',1.5);
axis tight; kgrid; kxlabel('time, $\mu$s');
kylabel('control amplitudes, MHz');
ktitle('band-limited optimal controls');
klegend({'cavity X','cavity Y','qubit X','qubit Y'},'Location','Best');

end

