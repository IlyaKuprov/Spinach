% 1H{14N} continuous-wave decoupling at 70 kHz MAS with the Ramsey
% shift broadening of the 1H line by the inhomogeneous 14N decoupler
% field: simulated reproduction of Figures 4 and 5a of Nehra, Agar-
% wal and Nishiyama, Solid State Nucl. Magn. Reson. 137 (2025) 102006
% for the amide 1H-14N pair of N-acetyl-L-alanine. 
%
% The residual dipolar splitting and the second-order quadrupolar
% shift of 14N are treated exactly using the numerical rotating fra-
% me transformation, with rotor-synchronous detection. The Ramsey 
% shift of the 1H resonance by the 14N decoupler field (Eq 14 of the
% manuscript) enters as a 1H chemical shift; B1 inhomogeneity of the
% probe is an explicit average over the B1 scaling factor sampled
% along the axis of a solenoid coil (+10% to -15% range used in the
% paper). Panel (a) shows the NH peak at each decoupling field (each
% spectrum is centred on its own maximum), panel (b) shows the drift
% and the broadening of the NH line with the decoupling field.
%
% Calculation time: about fifteen minutes on eleven workers
%
% ilya.kuprov@weizmann.ac.il

function ramsey_shift_14n_decoupling()

% Magnet field from the 1H frequency of 599.6 MHz
sys.magnet=2*pi*599.6e6/spin('1H');

% Larmor frequencies, Hz
nu0_h=sys.magnet*spin('1H')/(2*pi);
nu0_n=sys.magnet*spin('14N')/(2*pi);

% Amide 14N quadrupolar coupling, P_Q = 3.4 MHz
c_q=3.35e6; eta_q=0.30;

% 14N decoupling field strengths, Hz
nu1_n=(0:10)*4.2e3;

% B1 scaling factors along the axis of a solenoid coil
theta=linspace(-1,1,11)*acos(0.85/1.10);
b1_scale=1.10*cos(theta);

% Natural line width of the NH proton, Hz
nat_fwhm=200;

% Spin system: EFG z-axis normal to the amide plane, N-H (1.03 A) along EFG y-axis
sys.isotopes={'1H','14N'};
inter.coordinates={[0.00 0.00 0.00]; [0.00 1.03 0.00]};
inter.coupling.scalar=cell(2,2); inter.coupling.scalar{1,2}=65;
inter.coupling.matrix=cell(2,2);
inter.coupling.matrix{2,2}=eeqq2nqi(c_q,eta_q,1,[0 0 0]);

% Hilbert space formalism
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Hush the console
sys.output='hush';

% Experiment parameters
parameters.rate=70e3;
parameters.axis=[1 1 1];
parameters.max_rank=100;
parameters.grid='rep_2ang_200pts_sph';
parameters.rframes={{'14N',inf}};
parameters.spins={'1H','14N'};
parameters.npoints=1024;
parameters.zerofill=8192;
parameters.sweep=parameters.rate;
parameters.axis_units='ppm';
parameters.verbose=0;

% On-resonance 14N irradiation: carrier at the isotropic second-order shift
nu_q2=(3/40)*(c_q^2/nu0_n)*(1+eta_q^2/3);
parameters.offset=[7.96*nu0_h/1e6 nu_q2];

% Exponential apodisation matching the natural line width
k_apod=pi*nat_fwhm*parameters.npoints/parameters.rate;

% Loop over 14N decoupling field strengths
specs=zeros(parameters.zerofill,numel(nu1_n));
for n=1:numel(nu1_n)
    
    % Loop over B1 scaling factors
    for k=1:numel(b1_scale)

        % Ramsey shift of the 1H resonance by the 14N decoupler field
        nu1_h=b1_scale(k)*nu1_n(n)*spin('1H')/spin('14N');
        ramsey_hz=nu1_h^2*nu0_h/(nu0_h^2-nu0_n^2);

        % Observed 1H chemical shift includes the Ramsey shift
        inter.zeeman.scalar={7.96+1e6*ramsey_hz/nu0_h, 0.0};

        % Spinach housekeeping
        spin_system=create(sys,inter);
        spin_system=basis(spin_system,bas);

        % Initial and detection states
        parameters.rho0=state(spin_system,'L+','1H');
        parameters.coil=state(spin_system,'L+','1H');

        % 14N decoupling field for this B1 scaling factor
        parameters.rf_pwr=b1_scale(k)*nu1_n(n);

        % Simulation
        fid=singlerot(spin_system,@cw_acquire,parameters,'qnmr');

        % Apodisation, Fourier transform, and B1 averaging
        fid=apodisation(spin_system,fid,{{'exp',k_apod}});
        spec=real(fftshift(fft(fid,parameters.zerofill)));
        specs(:,n)=specs(:,n)+spec/numel(b1_scale);

    end

    % Report progress
    disp(['decoupling field ' num2str(nu1_n(n)/1e3) ' kHz done']);

end

% Chemical shift axis, referenced to the undecoupled NH resonance at 7.96 ppm
ax_pars=parameters; ax_pars.spins={'1H'}; ax_pars.offset=parameters.offset(1);
ppm=axis_1d(spin_system,ax_pars); ppm=ppm(:);
[~,idx]=max(specs(:,1)); ppm=ppm-ppm(idx)+7.96;

% NH peak intensities relative to the undecoupled spectrum
window=(ppm>6.5)&(ppm<9.5); heights=max(specs(window,:))/max(specs(window,1));

% Save the results
save('ramsey_shift_14n_decoupling.mat','ppm','specs','nu1_n','b1_scale','heights');

% Panel (a): NH peaks at each decoupling field, each centred on its own maximum
kfigure(); scale_figure([2.0 0.9]); subplot(1,2,1); hold on;
for n=1:numel(nu1_n)
    [~,idx]=max(specs(:,n).*window); mask=abs(ppm-ppm(idx))<1.5;
    plot(nu1_n(n)/1e3+1.4*(ppm(mask)-ppm(idx)),...
         specs(mask,n)/max(specs(window,1)),'b-');
end
plot(nu1_n/1e3,heights,'r-'); plot(nu1_n/1e3,heights,'ro');
kgrid; xlim([-2.1 44.1]); ylim([-0.05 1.4]);
kxlabel('$^{14}$N decoupling field strength, kHz'); 
kylabel('NH peak intensity, relative'); 
set(gca,'XTick',nu1_n/1e3); ktitle('(a)');

% Panel (b): NH line drift and broadening with the decoupling field
subplot(1,2,2); hold on; shown=1:2:numel(nu1_n);
for n=1:numel(shown)
    plot(ppm,specs(:,shown(n))/max(specs(window,shown(n)))+1.1*(n-1),'b-');
    text(6.6,1.1*(n-1)+0.55,['$\nu_{\rm 14N}$ = ' num2str(nu1_n(shown(n))/1e3) ' kHz'],...
         'Interpreter','latex','FontSize',11,'HorizontalAlignment','right');
end
plot([7.96 7.96],[-0.1 1.1*numel(shown)],'k--');
set(gca,'XDir','reverse','YTick',[]); xlim([6.5 10.0]); ylim([-0.1 1.1*numel(shown)]); 
kxlabel('$^{1}$H chemical shift, ppm'); ktitle('(b)'); kgrid;

% Export the figure
exportgraphics(gcf,'ramsey_shift_14n_decoupling.png','Resolution',300);
exportgraphics(gcf,'ramsey_shift_14n_decoupling.pdf','ContentType','vector');

end

% Rotor-synchronous acquisition under CW irradiation of 14N
function fid=cw_acquire(spin_system,parameters,H,~,~)

% Time step of the Hamiltonian stack
dt=1/(parameters.rate*numel(H));

% CW irradiation operator
H_rf=2*pi*parameters.rf_pwr*operator(spin_system,'Lx','14N');

% Rotor period propagator
P=eye(size(H{1}));
for n=1:numel(H)
    P=propagator(spin_system,H{n}+H_rf,dt)*P;
end

% Rotor-synchronous detection
fid=zeros(parameters.npoints,1); rho=parameters.rho0;
for n=1:parameters.npoints
    fid(n)=hdot(parameters.coil,rho); rho=P*rho*P';
end

end

