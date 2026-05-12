% Cavity-induced spin relaxation in the EPR Purcell regime.
% Coherent Jaynes-Cummings exchange is combined with rapid
% cavity damping in Liouville space, producing relaxation of
% the spin excitation by the NMR mechanism known as relaxation
% of the second kind.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function spin_cavity_purcell_effect()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'E','C3'};

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,[]);
spin_system=basis(spin_system,bas);

% Purcell parameters
coupling=2*pi*0.35e6;
detuning=2*pi*linspace(-8e6,8e6,301);
loss_rates=2*pi*[2e6 4e6 8e6 16e6];

% Spin detuning Hamiltonian
spin_ham=operator(spin_system,'Lz',1);

% Jaynes-Cummings exchange Hamiltonian
Hjc=coupling*(operator(spin_system,{'L+','A'},{1,2})+...
              operator(spin_system,{'L-','C'},{1,2}));

% Clean up numerical asymmetry
Hjc=(Hjc+Hjc')/2;

% Project the Hamiltonian coupling into the active doublet
spin_exc=state(spin_system,{'ZL2','BL1'},{1,2});
cav_exc=state(spin_system,{'ZL1','BL2'},{1,2});
jc_coupling=norm(cav_exc*Hjc*spin_exc,'fro');

% Validate the active matrix element
if abs(jc_coupling-coupling)>1e-10*coupling
    error('Jaynes-Cummings matrix element is inconsistent.');
end

% Build the cavity loss Lindblad dissipator in Liouville space
cav_ann=operator(spin_system,'A',2);
cav_pop=cav_ann'*cav_ann;
unit=speye(size(Hjc,1));
loss_gen=kron(conj(cav_ann),cav_ann)-...
         0.5*kron(unit,cav_pop)-...
         0.5*kron(cav_pop.',unit);

% Extract Purcell rates from Liouvillian slow modes
rates=zeros(numel(detuning),numel(loss_rates));
for n=1:numel(loss_rates)
    for k=1:numel(detuning)

        % Assemble the coherent Hamiltonian at this detuning
        H=detuning(k)*spin_ham+Hjc;
        H=(H+H')/2;

        % Assemble the Hamiltonian and dissipative Liouvillian
        L=-1i*(kron(unit,H)-kron(H.',unit))+loss_rates(n)*loss_gen;

        % Convert the slow spin-amplitude mode into a population rate
        decay_modes=eig(full(L));
        decay_modes=decay_modes(real(decay_modes)<-1e-3);
        rates(k,n)=-2*max(real(decay_modes));

    end
end

% Validate resonant second-kind relaxation
if rates(ceil(end/2),2)<=rates(1,2)
    error('Purcell rate is not resonantly enhanced.');
end

% Pick representative detunings for survival curves
time_axis=linspace(0,40e-6,250);
det_pick=2*pi*[0 2e6 6e6];
survival=zeros(numel(time_axis),numel(det_pick));
rho_vec=spin_exc(:);
spin_obs=spin_exc(:)';
loss_ref=2*pi*4e6;

% Simulate spin excitation survival under cavity damping
for n=1:numel(det_pick)

    % Assemble the coherent Hamiltonian at this detuning
    H=det_pick(n)*spin_ham+Hjc;
    H=(H+H')/2;

    % Assemble the Hamiltonian and dissipative Liouvillian
    L=-1i*(kron(unit,H)-kron(H.',unit))+loss_ref*loss_gen;

    % Propagate the density matrix in Liouville space
    for k=1:numel(time_axis)
        survival(k,n)=real(spin_obs*(expm(full(L)*time_axis(k))*rho_vec));
    end

end

% Validate the survival dynamics
if survival(end,1)>=survival(end,end)
    error('resonant Purcell relaxation is not faster than off-resonant relaxation.');
end

% Plot Liouvillian rates and survival dynamics
kfigure(); scale_figure([2.0 0.75]);
subplot(1,2,1); plot(detuning/(2*pi*1e6),rates/(2*pi*1e3),'LineWidth',1.5);
axis tight; kgrid; kxlabel('spin-cavity detuning, MHz');
kylabel('$\Gamma_P/2\pi$, kHz');
ktitle('Liouvillian Purcell rate');
klegend({'2 MHz','4 MHz','8 MHz','16 MHz'},'Location','NorthEast');
subplot(1,2,2); plot(1e6*time_axis,survival,'LineWidth',1.5);
axis tight; kgrid; kxlabel('time, $\mu$s');
kylabel('spin excitation survival');
ktitle('relaxation of the second kind');
klegend({'0 MHz','2 MHz','6 MHz'},'Location','NorthEast');

end
