% Simulation of nutation frequency dependent XiX DNP field profiles in 
% the steady state with electron-proton distance and electron Rabi 
% frequency ensembles.
%
% Calculation time: minutes
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de
% i.kuprov@soton.ac.uk

close all

nu=6.8;                                 % Electron nutation frequency, MHz
srt=0.051;                              % Repetition time, ms
for n=1
    xix_field_profile_b1_r(nu(n),srt(n))
end

nu=[9.6 13.5 18.5 25 36];               % Electron nutation frequency, MHz
srt=[0.051 0.102 0.153 0.153 0.306];    % Repetition time, ms
for n=1:numel(nu)
    xix_field_profile_b1_r(nu(n),srt(n))
end
legend('7 MHz','10 MHz','14 MHz','18 MHz','25 MHz','36 MHz')
set(gcf,'outerposition',[1000 650 710 510])

% Save results
savefig(gcf,'xix_nutation_ensemble_b1_r.fig');


function xix_field_profile_b1_r(nu,srt)

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

% Distance and B1 ensemble, Gauss-Legendre points
[r,wr]=gaussleg(3.5,20,3);      % Angstrom
if nu == 6.8
    [b1,wb1]=gaussleg(1e6,(nu+2)*1e6,5);  % Hz
else
    [b1,wb1]=gaussleg((nu-8)*1e6,(nu+2)*1e6,5);  % Hz
end

% Microwave resonance offsets, Hz
offsets=-52e6:-1e6:-64e6;

% Preallocate equilibrium DNP value array
dnp=zeros([numel(offsets) numel(r) numel(b1)],'like',1i);

% Over distances
for n=1:numel(r)

    % Cartesian coordinates
    inter.coordinates={[0.000 0.000 0.000];
                       [0.000 0.000 r(n) ]};
       
    % Relaxation rates, distance and orientation 
    % dependence provided using a function handle
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
    parameters.coil=state(spin_system,'Lz','1H');

    % Experiment parameters
    parameters.spins={'E','1H'};
    parameters.grid='rep_2ang_800pts_sph';
    parameters.pulse_dur=48e-9;              % Pulse duration, seconds
    parameters.nloops=36;                    % Number of XiX DNP blocks
    parameters.phase=pi;                     % Second pulse inverted phase
    parameters.addshift=-13e6;
    parameters.el_offs=offsets;

    % Calculate shot spacing
    parameters.shot_spacing=(srt*1e-3) - 2*parameters.nloops*parameters.pulse_dur;

    % Over B1 fields
    for k=1:numel(b1)

        % Set electron nutation frequency
        parameters.irr_powers=b1(k);

        % Run the steady state simulation
        dnp(:,n,k)=powder(spin_system,@xixdnp_steady,parameters,'esr');

    end

end

% Integrate over the B1 field distribution
dnp=sum(dnp.*reshape(wb1,[1 1 numel(wb1)]),3)/sum(wb1);

% Integrate over the distance distribution, r^2 is the radial part of the Jacobian
dnp=sum(dnp.*reshape(r.^2,[1 numel(r)]).*reshape(wr,[1 numel(wr)]),2)/sum((r.^2).*wr);

% Plotting 
figure(1); plot3(nu*ones(1,numel(offsets)),parameters.el_offs/1e6,-real(dnp),'o-','LineWidth',1.5); 
xlabel('\nu_1_S (MHz)');
ylabel('\Omega/2\pi (MHz)');
zlabel('\langle I_Z \rangle');
ax=gca;
ax.FontSize=14;
ax.LineWidth=1.2;
set(gca,'XMinorTick','on','YMinorTick','on');
grid on; hold on
view([-67 11]); xlim([0 40]); ylim([-65 -50]); zlim([0 1.2e-3]);

end

