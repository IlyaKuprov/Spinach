% Simulation of T1n dependent XiX DNP field profiles in the steady state 
% with electron-proton distance ensemble.
% 
% Calculation time: minutes
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de
% i.kuprov@soton.ac.uk

close all

xix_field_profile_ensemble_r_ori

T1n=[52 5 0.5 50e-3 5e-3];
Color={'#D95319' '#EDB120' '#77AC30' '#0072BD' '#7E2F8E'};

for j=1:numel(T1n)
    col=char(Color(j));
    xix_field_profile_ensemble_r(T1n(j),col)
end

legend('Ori. dep. {\itT}_1_n','52 s','5 s','0.5 s','50 ms','5 ms','location','southeast')

% Save figure
savefig(gcf,'xix_field_profile_ensemble_r_T1n.fig');


function xix_field_profile_ensemble_r_ori()

% Q-band magnet
sys.magnet=1.2142;

% Electron and proton
sys.isotopes={'E','1H'};

% Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]};
inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10]};

% Spin temperature
inter.temperature=80;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Propagator accuracy
sys.tols.prop_chop=1e-12;

% Algorithmic options
sys.disable={'hygiene'}';

% Distance ensemble, Gauss-Legendre points
[r,w]=gaussleg(3.5,20,3);

% Microwave resonance offsets, Hz
offsets=linspace(-100e6,100e6,201);

% Compute DNP at each distance
dnp=zeros([numel(offsets) numel(r)]);
for n=1:numel(r)

    % Cartesian coordinates
    inter.coordinates={[0.000 0.000 0.000];
                       [0.000 0.000 r(n) ]};
    
    % Relaxation rates, distance and ori. dep. R1n
    inter.relaxation={'t1_t2'};
    r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature,...
                                   2.00230,1e-3,52,r(n),bet); 
    inter.r1_rates={1e3 r1n_rate};
    inter.r2_rates={200e3 50e3};
    inter.rlx_keep='diagonal';
    inter.equilibrium='dibari';

    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Detect the proton
    parameters.coil=state(spin_system,'Lz',2);

    % Experiment parameters
    parameters.spins={'E','1H'};
    parameters.irr_powers=18e6;            % Electron nutation frequency [Hz]
    parameters.grid='rep_2ang_800pts_sph';
    parameters.pulse_dur=48e-9;              % Pulse duration, seconds
    parameters.nloops=36;                    % Number of XiX DNP blocks (power of 2)
    parameters.phase=pi;                     % Second pulse inverted phase
    parameters.addshift=-13e6;
    parameters.el_offs=offsets;

    % Calculate shot spacing
    parameters.shot_spacing=204e-6 - 2*parameters.nloops*parameters.pulse_dur;

    % Run the steady state simulation
    dnp(:,n)=powder(spin_system,@xixdnp_steady,parameters,'esr');

end

% Integrate over the distance distribution, r^2 is the Jacobian
dnp=sum(dnp.*reshape(r.^2,[1 numel(r)]).*reshape(w,[1 numel(w)]),2)/sum((r.^2).*w);

% Plotting 
figure(1); plot(parameters.el_offs/1e6,real(dnp),'color','#000000','LineWidth',1.5); 
xlabel('\Omega/2\pi (MHz)');
ylabel('\langle I_Z \rangle');
grid on; xlim tight; ylim([-1.5e-3 1.5e-3]); hold on

ax=gca;
ax.FontSize=14;
ax.LineWidth=1.2;
set(gca,'XMinorTick','on','YMinorTick','on');

end

function xix_field_profile_ensemble_r(T1n,col)

% Q-band magnet
sys.magnet=1.2142;

% Electron and proton
sys.isotopes={'E','1H'};

% Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]};
inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10]};

% Spin temperature
inter.temperature=80;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Propagator accuracy
sys.tols.prop_chop=1e-12;

% Algorithmic options
sys.disable={'hygiene'}';

% Distance ensemble, Gauss-Legendre points
[r,w]=gaussleg(3.5,20,3);

% Microwave resonance offsets, Hz
offsets=linspace(-100e6,100e6,201);

% Compute DNP at each distance
dnp=zeros([numel(offsets) numel(r)]);
for n=1:numel(r)

    % Cartesian coordinates
    inter.coordinates={[0.000 0.000 0.000];
                       [0.000 0.000 r(n) ]};
    
    % Relaxation rates, distance and ori. dep. R1n
    inter.relaxation={'t1_t2'};
    inter.r1_rates={1e3 1/T1n};
    inter.r2_rates={200e3 50e3};
    inter.rlx_keep='diagonal';
    inter.equilibrium='dibari';

    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Detect the proton
    parameters.coil=state(spin_system,'Lz',2);

    % Experiment parameters
    parameters.spins={'E','1H'};
    parameters.irr_powers=18e6;            % Electron nutation frequency [Hz]
    parameters.grid='rep_2ang_800pts_sph';
    parameters.pulse_dur=48e-9;              % Pulse duration, seconds
    parameters.nloops=36;                    % Number of XiX DNP blocks (power of 2)
    parameters.phase=pi;                     % Second pulse inverted phase
    parameters.addshift=-13e6;
    parameters.el_offs=offsets;

    % Calculate shot spacing
    parameters.shot_spacing=204e-6 - 2*parameters.nloops*parameters.pulse_dur;

    % Run the steady state simulation
    dnp(:,n)=powder(spin_system,@xixdnp_steady,parameters,'esr');

end

% Integrate over the distance distribution, r^2 is the Jacobian
dnp=sum(dnp.*reshape(r.^2,[1 numel(r)]).*reshape(w,[1 numel(w)]),2)/sum((r.^2).*w);

% Plotting 
figure(1); plot(parameters.el_offs/1e6,real(dnp),'color',col,'LineWidth',1.5); 
xlabel('\Omega/2\pi (MHz)');
ylabel('\langle I_Z \rangle');
grid on; xlim tight; ylim([-1.5e-3 1.5e-3]); hold on

ax=gca;
ax.FontSize=14;
ax.LineWidth=1.2;
set(gca,'XMinorTick','on','YMinorTick','on');

end

