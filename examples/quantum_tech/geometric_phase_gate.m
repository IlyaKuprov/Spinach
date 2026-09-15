% Geometric phase gate between two trapped-ion qubits driven by a
% state-dependent optical dipole force, as demonstrated by Leibfried
% et al. (Nature 422, 412 (2003)) with two 9Be+ ions. Here a three-
% ion chain in the same trap is simulated: a global force beam pair
% detuned by delta from the stretch mode pushes the outer ions in op-
% posite directions, while the middle ion has zero amplitude in the
% stretch mode eigenvector and so does not couple to that mode; the
% outer ions acquire a state-dependent geometric phase and the middle
% ion is a spectator. All three ions also couple to the centre-of-mass
% and Egyptian modes of the chain (frequencies from James, Appl. Phys.
% B 66, 181 (1998)), which are far off resonance with the drive but are
% kept in the simulation with their proper force constants. The ion qubits are
% treated as spin-1/2 particles, the motional modes as bosonic modes,
% and the Hamiltonian is written in the frame rotating with the force
% drive, where the spin-dependent force is a static longitudinal cou-
% pling and the mode frequencies are their detunings from the drive,
% supplied as parameters.mode_offset of the device context.
%
% The basis set is IK-SBS: correlations are traced separately on the
% boson-boson, spin-boson, and spin-spin coupling graphs. Pure spin-
% spin correlations up to order two are kept inside the spin-boson
% subgraphs, which is what the entangling phase requires.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function geometric_phase_gate()

% Trap parameters, Leibfried et al. (2003)
stretch_frq=6.1e6;         % stretch mode frequency, Hz
delta=26e3;                % force drive detuning from the stretch mode, Hz

% Other normal modes of a three-ion chain, James (1998)
com_frq=stretch_frq/sqrt(3);
egypt_frq=sqrt(29/5)*com_frq;

% Force drive frequency
drive_frq=stretch_frq+delta;

% Stretch mode force constant giving a pi/2 geometric phase at T=1/delta
kappa=delta/4;

% Lamb-Dicke scaling of the force constants of the other modes
kappa_com=kappa*sqrt(2/3)*sqrt(stretch_frq/com_frq);
kappa_egypt=kappa*sqrt(2/6)*sqrt(stretch_frq/egypt_frq);

% Ion qubits and the three axial modes of the chain
sys.magnet=0;
sys.isotopes={'E','E','E','V6','V3','V3'};

% Laboratory frame mode frequencies
inter.modes.frqs={[] [] [] stretch_frq com_frq egypt_frq};

% Spin-dependent forces, Spinach convention Lz(a+a')/sqrt(2)
inter.modes.longitudinal=cell(6,6);
inter.modes.longitudinal{1,4}=+2*sqrt(2)*kappa;
inter.modes.longitudinal{3,4}=-2*sqrt(2)*kappa;
inter.modes.longitudinal{1,5}=+2*sqrt(2)*kappa_com;
inter.modes.longitudinal{2,5}=+2*sqrt(2)*kappa_com;
inter.modes.longitudinal{3,5}=+2*sqrt(2)*kappa_com;
inter.modes.longitudinal{1,6}=+2*sqrt(2)*kappa_egypt;
inter.modes.longitudinal{2,6}=-4*sqrt(2)*kappa_egypt;
inter.modes.longitudinal{3,6}=+2*sqrt(2)*kappa_egypt;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-SBS';
bas.connectivity='full_tensors';
bas.inter_level=[2 3 2];

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.rho0=state(spin_system,{'ZL1','ZL1','ZL1','BL1','BL1','BL1'},{1,2,3,4,5,6});
parameters.mode_offset=[drive_frq drive_frq drive_frq];
parameters.npoints=1561;
parameters.gate_time=1/delta;

% Run the gate through the device context in the frame of the force drive
answer=device(spin_system,@gate_sequence,parameters,'spin-phonon');

% Validate the spectator ion and the closure of the motional loop
if abs(answer.sigma_x(2,end)-answer.sigma_x(2,1))>1e-2
    error('spectator ion coherence is not preserved.');
end
if abs(answer.stretch_pop(end))>1e-2
    error('stretch mode did not return to its ground state.');
end

% Validate the entangling phase: parity contrast of the outer ions
if abs((max(answer.parity)-min(answer.parity))/2-1)>2e-2
    error('parity contrast is not consistent with a pi/2 geometric phase.');
end

% Plot the results
time_axis=1e6*linspace(0,parameters.gate_time,parameters.npoints);
kfigure(); scale_figure([2.0 0.75]);
subplot(1,2,1); plot(time_axis,answer.sigma_x,'LineWidth',1.5);
hold on; plot(time_axis,answer.stretch_pop,'LineWidth',1.5);
axis tight; kgrid; kxlabel('time, $\mu$s');
klegend({'$\sigma_x$, ion 1','$\sigma_x$, ion 2','$\sigma_x$, ion 3','$\langle a^{+}a\rangle$, stretch'},'Location','best');
subplot(1,2,2); plot(answer.phases/pi,answer.parity,'LineWidth',1.5);
axis tight; kgrid; kxlabel('analysis pulse phase, $\pi$ rad');
kylabel('parity of ions 1 and 3');

end

% State-dependent displacement between two global pi/2 pulses, then a parity scan
function answer=gate_sequence(spin_system,parameters,H,R,K)

% Compose the Liouvillian
L=H+1i*R+1i*K;

% Global rotation generators, observables, and the unit state for normalisation
Lx=operator(spin_system,'Lx','E'); Ly=operator(spin_system,'Ly','E');
coil_x=[state(spin_system,{'Lx'},{1}) state(spin_system,{'Lx'},{2}) state(spin_system,{'Lx'},{3})];
coil_n=state(spin_system,{'N'},{4});
coil_p=state(spin_system,{'Lz','Lz'},{1,3}); unit=state(spin_system,{'E'},{1});

% First pi/2 pulse puts every ion along +x
rho=step(spin_system,Ly,parameters.rho0,pi/2);

% State-dependent displacement over one closed loop of the stretch mode
dt=parameters.gate_time/(parameters.npoints-1);
traj=evolution(spin_system,L,[],rho,dt,parameters.npoints-1,'trajectory');
answer.sigma_x=real(2*(coil_x'*traj)./(ones(3,1)*(unit'*traj)));
answer.stretch_pop=real((coil_n'*traj)./(unit'*traj));

% Second pi/2 pulse turns the phase-gated state into a GHZ-type state of the outer ions
rho=step(spin_system,Ly,traj(:,end),pi/2);

% Parity of the outer ions after an analysis pi/2 pulse of variable phase
answer.phases=linspace(0,2*pi,101); answer.parity=zeros(size(answer.phases));
for n=1:numel(answer.phases)
    rho_an=step(spin_system,cos(answer.phases(n))*Lx+sin(answer.phases(n))*Ly,rho,pi/2);
    answer.parity(n)=real(4*(coil_p'*rho_an)/(unit'*rho_an));
end

end

