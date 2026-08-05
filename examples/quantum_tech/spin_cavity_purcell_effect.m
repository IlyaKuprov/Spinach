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
bas.formalism='zeeman-liouv';
bas.approximation='none';

% Purcell parameters
coupling=0.35e6;
detuning=2*pi*linspace(-8e6,8e6,301);
loss_rates=2*pi*[2e6 4e6 8e6 16e6];

% Preallocate rate array
rates=zeros(numel(detuning),numel(loss_rates));

% Loop over cavity loss rates
for n=1:numel(loss_rates)

    % Generators of the resonant damped spin-cavity device
    [spin_system,Hjc,R]=purcell_device(sys,bas,coupling,loss_rates(n));

    % Spin detuning superoperator
    spin_ham=operator(spin_system,'Lz',1);

    % Extract Purcell rates from Liouvillian slow modes
    for k=1:numel(detuning)

        % Assemble the dissipative Liouvillian at this detuning
        L=-1i*(Hjc+detuning(k)*spin_ham)+R;

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

% Validate the resonant rate against the exact two-level Purcell rate
g_rad=2*pi*coupling; kappa_ref=loss_rates(2);
purcell=kappa_ref/2-sqrt(kappa_ref^2/4-4*g_rad^2);
if abs(rates(ceil(end/2),2)-purcell)>1e-4*purcell
    error('resonant Purcell rate does not match the exact expression.');
end

% Rebuild the generators at the reference loss rate
[spin_system,Hjc,R]=purcell_device(sys,bas,coupling,loss_rates(2));
spin_ham=operator(spin_system,'Lz',1);

% Pick representative detunings for survival curves
time_axis=linspace(0,40e-6,250);
det_pick=2*pi*[0 2e6 6e6];
survival=zeros(numel(time_axis),numel(det_pick));
rho_vec=state(spin_system,{'ZL2','BL1'},{1,2});
spin_obs=rho_vec';

% Simulate spin excitation survival under cavity damping
for n=1:numel(det_pick)

    % Assemble the dissipative Liouvillian at this detuning
    L=-1i*(Hjc+det_pick(n)*spin_ham)+R;

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
klegend({'$\kappa/2\pi$=2 MHz','$\kappa/2\pi$=4 MHz',...
         '$\kappa/2\pi$=8 MHz','$\kappa/2\pi$=16 MHz'},'Location','Best');
subplot(1,2,2); plot(1e6*time_axis,survival,'LineWidth',1.5);
axis tight; kgrid; kxlabel('time, $\mu$s');
kylabel('spin excitation survival');
ktitle('Purcell decay curves');
klegend({'0 MHz','2 MHz','6 MHz'},'Location','Best');

end

% Generators of a resonant damped spin-cavity device
function [spin_system,Hjc,R]=purcell_device(sys,bas,coupling,loss_rate)
inter.modes.frqs={[] 0};
inter.modes.linewidths={[] loss_rate/(2*pi)};
inter.modes.exchange=cell(2,2);
inter.modes.exchange{1,2}=coupling;
inter.temperature=0;
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
Hjc=hamiltonian(assume(spin_system,'cavity'));
R=relaxation(spin_system);
end

