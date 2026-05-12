% Cavity-induced spin relaxation in the EPR Purcell regime.
% The Lorentzian Purcell rate follows Bienfait et al., Nature
% 531, 74-77 (2016), and is plotted together with the corres-
% ponding survival dynamics.
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
kappas=2*pi*[0.5e6 1.0e6 2.0e6 4.0e6];

% Jaynes-Cummings Hamiltonian in the one-excitation manifold
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

% Compute the Lorentzian Purcell rates
rates=zeros(numel(detuning),numel(kappas));
for n=1:numel(kappas)
    rates(:,n)=kappas(n)*jc_coupling^2./...
               (detuning.^2+(kappas(n)/2)^2);
end

% Validate resonant enhancement of the Purcell rate
if rates(ceil(end/2),2)<=rates(1,2)
    error('Purcell rate is not resonantly enhanced.');
end

% Pick representative detunings for survival curves
time_axis=linspace(0,40e-6,250);
det_pick=2*pi*[0 2e6 6e6];
survival=zeros(numel(time_axis),numel(det_pick));
for n=1:numel(det_pick)
    rate=kappas(2)*jc_coupling^2/(det_pick(n)^2+(kappas(2)/2)^2);
    survival(:,n)=exp(-rate*time_axis);
end

% Plot Purcell rate and survival dynamics
kfigure(); scale_figure([2.0 0.75]);
subplot(1,2,1); plot(detuning/(2*pi*1e6),rates/(2*pi*1e3),'LineWidth',1.5);
axis tight; kgrid; kxlabel('spin-cavity detuning, MHz');
kylabel('$\Gamma_P/2\pi$, kHz');
ktitle('Purcell rate');
klegend({'0.5 MHz','1 MHz','2 MHz','4 MHz'},'Location','NorthEast');
subplot(1,2,2); plot(1e6*time_axis,survival,'LineWidth',1.5);
axis tight; kgrid; kxlabel('time, $\mu$s');
kylabel('spin excitation survival');
ktitle('Purcell relaxation');
klegend({'0 MHz','2 MHz','6 MHz'},'Location','NorthEast');

end
