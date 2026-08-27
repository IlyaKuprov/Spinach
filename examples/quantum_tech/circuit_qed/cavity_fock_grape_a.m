% GRAPE preparation of a cavity Fock state through a dispersively
% coupled qubit, using piecewise-constant drives on both the cavity
% and the qubit. A linear drive alone cannot make a Fock state out
% of the vacuum of a harmonic mode; the qubit conditions the cavity
% phase through the dispersive shift and thereby provides the requi-
% red nonlinearity. The optimisation is run at two Fock space trun-
% cations; the pulse optimised in the smaller space underperforms
% when it is re-evaluated in the larger one - optimal control solu-
% tions must be converged with respect to the Fock space truncation.
% Model and parameters from the bosonic GRAPE example of the para-
% qeet package.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function cavity_fock_grape_a()

% Cavity truncation levels to compare
trunc_labels={'C3','C4'};

% Preallocate pulse and fidelity storage
pulses=cell(1,2); fids=zeros(1,3);

% Loop over the cavity truncations
for m=1:2

    % Magnet field
    sys.magnet=0;

    % Truncated cavity mode and a qubit
    sys.isotopes={trunc_labels{m},'E'};

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

    % Define control parameters
    control.isotopes={trunc_labels{m},'E'};
    control.channels=[1;1;2;2];
    control.drifts={{H}};
    control.operators=ops;
    control.rho_init={rho_init};
    control.rho_targ={rho_targ};
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

    % Gaussian initial guess on the in-phase channels
    env=exp(-(33e-9*((1:40)-0.5)-0.66e-6).^2/(2*(1.32e-6/8)^2));
    guess=[0.7*env; 0*env; 0.7*env; 0*env];

    % Run the optimisation, get normalised pulse
    pulses{m}=fmaxnewton(spin_system,@grape_xy,guess);

    % Recompute the transfer fidelity by direct propagation
    fids(m)=pulse_fidelity(spin_system,H,ops,pulses{m},rho_init,rho_targ);

end

% Re-evaluate the small truncation pulse in the larger space
fids(3)=pulse_fidelity(spin_system,H,ops,pulses{1},rho_init,rho_targ);

% Validate the optimisations
if (fids(1)<0.95)||(fids(2)<0.95)
    error('GRAPE optimisation did not converge.');
end

% Report the Fock truncation convergence lesson
disp(['Fidelity lost by transplanting the N=3 pulse into N=4: ' ...
      num2str(fids(2)-fids(3))]);

% Plot the fidelity comparison
kfigure(); bar(fids); kgrid;
set(gca,'XTickLabel',{'opt & eval N=3','opt & eval N=4','opt N=3, eval N=4'});
kylabel('Fock state transfer fidelity');
ktitle('Fock space truncation convergence');

end

% Transfer fidelity of a piecewise-constant pulse
function fid=pulse_fidelity(spin_system,H,ops,pulse,rho_init,rho_targ)
rho=rho_init;
for n=1:size(pulse,2)
    slice_ham=H+spin_system.control.pwr_levels*...
                (pulse(1,n)*ops{1}+pulse(2,n)*ops{2}+...
                 pulse(3,n)*ops{3}+pulse(4,n)*ops{4});
    P=propagator(spin_system,slice_ham,spin_system.control.pulse_dt(n));
    rho=P*rho*P';
end
fid=real(trace(rho_targ'*rho));
end

