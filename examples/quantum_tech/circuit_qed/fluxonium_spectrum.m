% Energy spectrum of a fluxonium qubit as a function of the exter-
% nal flux, reproducing Fig. 2.3(b,c,d) of Y. Lu, Optimal Control
% and Coherence Engineering for Superconducting Qubits, PhD thesis,
% Northwestern University, 2026. The Hamiltonian is built in the
% truncated harmonic oscillator basis by fluxonium.m for EJ/h=5.0
% GHz, EC/h=1.0 GHz, and EL/h=1.0 GHz; the three lowest levels are
% plotted against the external flux, and the potential landscape
% with the two lowest probability densities is shown at 0.4 and
% 0.5 flux quanta. The harmonic limit, the flux periodicity, the
% half-flux sweet spot, the Hermiticity of the Hamiltonian, and
% the convergence with respect to the basis size are validated.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function fluxonium_spectrum()

% Circuit energies, Hz
ec=1.0e9; ej=5.0e9; el=1.0e9;

% Flux sweep and oscillator basis sizes
flux_pts=linspace(0,1,201); nlev=[30 60];

% Three lowest levels in GHz for both basis sizes
levels=zeros(3,numel(flux_pts),2);
for m=1:2
    for k=1:numel(flux_pts)
        H=fluxonium(ec,ej,el,2*pi*flux_pts(k),nlev(m));
        energies=eig(H); levels(:,k,m)=energies(1:3)/(2*pi*1e9);
    end
end

% Hermiticity check on the last Hamiltonian
if norm(H-H',1)>1e-12*norm(H,1)
    error('fluxonium Hamiltonian is not Hermitian.');
end

% Qubit frequency from the bigger basis
f01=levels(2,:,2)-levels(1,:,2);

% Report the qubit frequency at zero and half flux
disp(['f01 at zero flux: ' num2str(f01(1),'%.6f') ' GHz']);
disp(['f01 at half flux: ' num2str(f01(101),'%.6f') ' GHz']);

% Convergence of the levels and of their spacings
lev_conv=max(abs(levels(:,:,1)-levels(:,:,2)),[],'all')/max(abs(levels(:,:,2)),[],'all');
spacings=levels(2:3,:,:)-levels([1 1],:,:);
spc_conv=max(abs(spacings(:,:,1)-spacings(:,:,2))./spacings(:,:,2),[],'all');
disp(['level convergence, 30 vs 60 states: ' num2str(lev_conv,'%.2e')]);
disp(['spacing convergence, 30 vs 60 states: ' num2str(spc_conv,'%.2e')]);
if (lev_conv>1e-6)||(spc_conv>1e-6)
    error('levels are not converged with respect to the basis size.');
end

% Harmonic limit: tiny Josephson energy gives the plasma frequency
energies=eig(fluxonium(ec,1,el,0,nlev(2)));
harm_err=abs((energies(2)-energies(1))/(2*pi)-sqrt(8*ec*el))/sqrt(8*ec*el);
disp(['harmonic limit deviation: ' num2str(harm_err,'%.2e')]);
if harm_err>1e-6
    error('harmonic limit does not reproduce the plasma frequency.');
end

% Sweet spot at half flux and periodicity in the flux quantum
[~,min_idx]=min(f01);
if min_idx~=101
    error('qubit frequency is not minimal at half flux.');
end
if max(abs(levels(:,1,2)-levels(:,end,2)))>1e-6*max(abs(levels(:,1,2)))
    error('spectrum is not periodic in the flux quantum.');
end

% Spectrum as a function of the external flux
kfigure(); scale_figure([2.4 0.75]);
subplot(1,3,1); plot(flux_pts,levels(:,:,2)','LineWidth',1.5);
axis tight; kgrid; kxlabel('$\Phi_{ext}/\Phi_0$'); kylabel('energy/h, GHz');
ktitle('fluxonium spectrum');

% Phase grid and normalised Hermite functions of the oscillator basis
phi_grid=linspace(-5,5,501); phi_zpf=(2*ec/el)^(1/4); xi=phi_grid/phi_zpf;
herm=zeros(nlev(2),numel(phi_grid));
herm(1,:)=exp(-xi.^2/2)/(pi^(1/4)*sqrt(phi_zpf));
herm(2,:)=sqrt(2)*xi.*herm(1,:);
for k=2:(nlev(2)-1)
    herm(k+1,:)=sqrt(2/k)*xi.*herm(k,:)-sqrt((k-1)/k)*herm(k-1,:);
end

% Potential and the two lowest probability densities at 0.4 and 0.5 flux quanta
for m=1:2
    flux=0.3+0.1*m; phi_e=2*pi*flux;
    H=fluxonium(ec,ej,el,phi_e,nlev(2)); [V,D]=eig(H);
    potential=(-ej*cos(phi_grid-phi_e)+(el/2)*phi_grid.^2)/1e9;
    densities=abs(herm'*V(:,1:2)).^2; energies=diag(D)/(2*pi*1e9);
    subplot(1,3,m+1); plot(phi_grid,potential,'Color',[0.5 0.5 0.5],'LineWidth',1.5);
    hold on; plot(phi_grid,energies(1)+densities(:,1),'LineWidth',1.5);
    plot(phi_grid,energies(2)+densities(:,2),'LineWidth',1.5); ylim([-5 5]);
    kgrid; kxlabel('$\varphi$, rad'); kylabel('energy/h, GHz');
    ktitle(['$\Phi_{ext}=' num2str(flux) '\Phi_0$']);
    klegend({'potential','$|\psi_0|^2$','$|\psi_1|^2$'},'Location','North');
end

end


