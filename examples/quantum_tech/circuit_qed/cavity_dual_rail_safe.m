% Flux-noise dephasing rates of two dual-rail qubits whose flux-tuna-
% ble transmon-coupled rails share one transmon ancilla, and their
% suppression by a Stark-assisted flux-noise evasion (SAFE) drive on
% the transmon, Sec. 4.4.2 and Fig. 4.4(c) of Yunwei Lu's PhD thesis
% (Northwestern University, 2026). The logical dephasing rate of each
% dual-rail qubit is set by the flux sensitivity of the dressed sin-
% gle-photon transition frequency of its transmon-coupled rail, Eqs.
% (4.93) and (4.95); the rates are computed from the eigenvalues of
% the rotating frame Hamiltonian of Eq. (4.16) at a fixed drive amp-
% litude as functions of the transmon-drive detuning. Both rates go
% through a minimum in the same detuning window, so that one drive
% protects both qubits.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function cavity_dual_rail_safe()

% Transmon anharmonicity, Hz
anharm=-200e6;

% Exchange couplings and detunings of the two transmon-coupled rails, Hz
g_bc=[50e6 60e6]; delta_bc=[1.414e9 2.0e9];

% Dispersive shifts and the sensitivities of the rail terms to the transmon frequency
chi=2*(g_bc./delta_bc).^2*anharm; sens_c=(g_bc./delta_bc).^2; sens_x=-4*anharm*g_bc.^2./delta_bc.^3;

% Transmon frequency sensitivity to flux, rad/s per flux quantum, and the noise amplitude, flux quanta
dwb_dphi=2*pi*6e9; noise_amp=1e-5;

% Drive amplitude, rad/s, and the transmon-drive detunings to scan, Hz
omega0=2*pi*10e6; detunings=-(20:0.5:50)*1e6;

% Three-level transmon and the two transmon-coupled rails, all in their own frames
sys.magnet=0; sys.isotopes={'T3','C3','C3'};
inter.modes.frqs={0 0 0};
inter.modes.anharms={anharm [] []};
inter.modes.kerr=cell(3,3); inter.modes.kerr{1,2}=chi(1); inter.modes.kerr{1,3}=chi(2);
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Hamiltonian without the detuning term, detuning, drive, and flux noise operators
H_hilb=hamiltonian(assume(spin_system,'labframe'));
num_b=operator(spin_system,'N',1);
drive_op=operator(spin_system,'C',1)+operator(spin_system,'A',1);
noise_op=num_b+sens_c(1)*operator(spin_system,'N',2)+sens_c(2)*operator(spin_system,'N',3)+...
       sens_x(1)*operator(spin_system,{'N','N'},{1,2})+sens_x(2)*operator(spin_system,{'N','N'},{1,3});

% Positions of the bare |g,0,0>, |g,1,0>, and |g,0,1> states
idx=[find(diag(state(spin_system,{'BL1','BL1','BL1'},{1,2,3}))>0.5) ...
     find(diag(state(spin_system,{'BL1','BL2','BL1'},{1,2,3}))>0.5) ...
     find(diag(state(spin_system,{'BL1','BL1','BL2'},{1,2,3}))>0.5)];

% Flux noise dephasing rates of the two rails under the drive, Eq. (4.93) with |ln(omega_ir*t)|=4
rates=zeros(2,numel(detunings)); dw=2*pi*1e3; shifts=[-dw dw];
for n=1:numel(detunings)
    energies=zeros(2,3);
    for m=1:2
        [vecs,vals]=eig(full(H_hilb+2*pi*detunings(n)*num_b+omega0*drive_op+shifts(m)*noise_op));
        for j=1:3
            [~,col]=max(abs(vecs(idx(j),:))); energies(m,j)=vals(col,col);
        end
    end
    suscept=(energies(2,2:3)-energies(2,1)-energies(1,2:3)+energies(1,1))/(2*dw);
    rates(:,n)=noise_amp*dwb_dphi*sqrt(2*4)*abs(suscept)';
end

% Undriven rates from the bare dispersive sensitivities, Eq. (4.95)
rates_off=noise_amp*dwb_dphi*sqrt(2*4)*sens_c;

% Locate the minima and the common operating point
[~,min_idx]=min(rates,[],2); common=round(mean(min_idx));
disp(['rail 1 minimum at ' num2str(-detunings(min_idx(1))/1e6,'%.1f') ' MHz, rail 2 minimum at ' ...
      num2str(-detunings(min_idx(2))/1e6,'%.1f') ' MHz, common point ' num2str(-detunings(common)/1e6,'%.1f') ' MHz']);
disp(['suppression factors at the common point: ' num2str(rates_off(1)/rates(1,common),'%.1f') ...
      ' and ' num2str(rates_off(2)/rates(2,common),'%.1f')]);

% Validate the sweet spots
if any(min_idx==1)||any(min_idx==numel(detunings))
    error('dephasing rate minima are outside the detuning window.');
end
if abs(detunings(min_idx(1))-detunings(min_idx(2)))>5e6
    error('the two rails do not share a sweet spot window.');
end
if any(rates_off'./rates(:,common)<5)
    error('SAFE drive does not suppress both dephasing rates fivefold.');
end

% Plot the dephasing rates, Fig. 4.4(c)
kfigure(); semilogy(-detunings/1e6,1e-9*rates',-detunings([1 end])/1e6,1e-9*[1; 1]*rates_off,'--'); kgrid;
kxlabel('$-\Delta_{bd}/2\pi$, MHz'); kylabel('dephasing rate, 1/ns');
klegend({'rail 1, driven','rail 2, driven','rail 1, undriven','rail 2, undriven'},'Location','northeast');
ktitle('dual-rail qubits under SAFE');

end

