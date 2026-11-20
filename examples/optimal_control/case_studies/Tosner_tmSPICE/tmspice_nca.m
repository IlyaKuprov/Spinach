% Reduced tmSPICE N-to-CA transfer replay after Tosner et al. (2018).
% Uses the four-spin model and two sampled trajectories from the 3.2 mm
% probe map. The tabulated spatial weights already include carbon receive
% magnitude. This short replay is not the full powder or pulse design.
%
% Talos <talos@spindynamics.org>

function tmspice_nca()

% Specify the four-spin model at 850 MHz proton frequency
sys.magnet=2*pi*850e6/spin('1H');
sys.isotopes={'15N','13C','13C','13C'};
sys.parallel={'processes',2};
sys.output='hush';
iso=[0 0 118.4 -16.1];
aniso=[99 -20 -76 -20];
asym=[0.19 0.43 0.90 0.43];
angles=[ 103.01350 -141.57466  65.13328;
         -81.06377   37.80482  37.44094;
        -150.38506   89.66823 -24.46362;
         -81.06377   37.80482  37.44094];
for n=1:4

    % Convert the shift parameters into principal values in ppm
    inter.zeeman.eigs{n}=iso(n)+aniso(n)*[-(1+asym(n))/2 (asym(n)-1)/2 1];
    inter.zeeman.euler{n}=deg2rad(angles(n,:));

end

% Use the six dipolar tensors and three scalar couplings
pairs=[1 2;1 3;1 4;2 3;2 4;3 4];
dip=[1019.17737 224.75205 207.87002 -2153.96075 -2119.03837 -476.30486];
angles=[63.76109 113.83214;92.27945 88.65551;75.88146 149.71407;...
        61.21563 -115.71998;88.50250 -0.56171;105.97922 29.37194];
inter.coupling.matrix=cell(4);
for n=1:6

    % Orient the dipolar axis in the molecular frame
    D=euler2dcm(deg2rad(angles(n,2)),deg2rad(angles(n,1)),0);
    inter.coupling.matrix{pairs(n,1),pairs(n,2)}=dip(n)*D*diag([-1 -1 2])*D';

end
inter.coupling.scalar=cell(4);
inter.coupling.scalar{1,2}=-11;
inter.coupling.scalar{2,3}=55;
inter.coupling.scalar{2,4}=35;

% Spinach housekeeping
bas.formalism='zeeman-hilb';
bas.approximation={'none'};
spin_system=basis(create(sys,inter),bas);
spin_system=assume(spin_system,'labframe');

% Use two orientations, two initial phases, and six midpoints per revolution
parameters.grid=struct('alphas',[0 0],'betas',[0.45 1.15],...
                       'gammas',[0 1.7],'weights',[0.4 0.6]);
parameters.rframes={{'15N',1},{'13C',1}};
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.n_ticks=6;
parameters.n_phases=2;
parameters.n_slices=12;
[drifts,weights]=mas_drifts(spin_system,parameters);

% Specify the two active RF channels and normalised transfer operators
control.isotopes={'15N','13C'};
control.channels=[1;1;2;2];
control.operators={operator(spin_system,'Lx','15N'),operator(spin_system,'Ly','15N'),...
                   operator(spin_system,'Lx','13C'),operator(spin_system,'Ly','13C')};
rho=state(spin_system,'Lx',1);
target=coil_state(spin_system,'Lx',2);
control.rho_init={rho/norm(rho,'fro')};
control.rho_targ={target/norm(target,'fro')};
control.pwr_levels=2*pi*30e3;
control.pulse_dt=ones(1,12)/(6*16500);
control.drifts=drifts;
control.drift_weights=weights;
control.method='lbfgs';
control.max_iter=2;
control.plotting={};
control.penalties={'none'};
control.p_weights=0;
control.l_bound=-1;
control.u_bound=1;

% Seed a reproducible waveform and verify the unity-map limit
rng(2018);
wave=0.25*randn(4,12);
[~,bare]=grape_xy(wave,optimcon(spin_system,control));
control.distortion={@(w)rotor_rf(w,ones(2,12))};
[~,unity]=grape_xy(wave,optimcon(spin_system,control));
assert(abs(bare(1)-unity(1))<1e-12);

% Use two six-sample probe trajectories and their receive-inclusive weights
maps={[0.442773-0.00319435i 0.564436-0.196016i 0.648825-0.201133i 0.655091+0.000655143i 0.536485+0.189028i 0.393987+0.178788i;...
       0.44231-0.0031699i 0.56383-0.195844i 0.648218-0.200969i 0.654553+0.000649118i 0.536091+0.188894i 0.393699+0.178681i];...
      [0.560875+0.0224739i 0.705408-0.128776i 0.792774-0.143677i 0.813989+0.0206696i 0.697284+0.215611i 0.525179+0.213355i;...
       0.560318+0.0224477i 0.704678-0.128774i 0.79207-0.143651i 0.81335+0.0206307i 0.696781+0.215473i 0.5248+0.213234i]};
receive=[0.0580128 0.0720874];

% Replicate the physical drift outside the rotor stack and link the maps
control.drifts=repelem(drifts,2);
control.drift_weights=kron(weights,receive(:));
control.distortion=cell(numel(control.drifts),1);
sampled=cell(size(control.distortion));
for n=1:numel(drifts)

    % Interpolate the periodic samples at the drift tick midpoints
    tick=mod((0:11)+3*mod(n-1,2),6)+0.5;
    for k=1:2

        % Reuse the sampled gain array for energy evaluation
        gain=interp1(0:6,[maps{k} maps{k}(:,1)].',tick,'linear').';
        idx=2*(n-1)+k;
        sampled{idx}=gain;
        control.distortion{idx}=@(w)rotor_rf(w,gain);

    end

end
control.ens_corrs={'dist_drift'};
control.off_ops={operator(spin_system,'Lz','15N'),operator(spin_system,'Lz','13C')};
control.offsets={[-500 500],[-600 600]};
spin_system=optimcon(spin_system,control);

% Verify a directional derivative on the linked, receive-weighted ensemble
[~,before,grad]=grape_xy(wave,spin_system);
direction=randn(size(wave));
direction=direction/norm(direction,'fro');
delta=1e-5;
[~,plus]=grape_xy(wave+delta*direction,spin_system);
[~,minus]=grape_xy(wave-delta*direction,spin_system);
slope=sum(grad(:,:,1).*direction,'all');
fd=(plus(1)-minus(1))/(2*delta);
assert(abs(slope-fd)<1e-6);

% Optimise briefly and evaluate energy on the same spatial ensemble
wave=fmaxnewton(spin_system,@grape_xy,wave);
[~,after]=grape_xy(wave,spin_system);
energy=zeros(size(sampled));
for n=1:numel(energy)

    % Integrate the squared angular-frequency drive over the pulse
    rf=control.pwr_levels*rotor_rf(wave,sampled{n});
    energy(n)=sum(sum(rf.^2,1).*control.pulse_dt);

end
avg_energy=dot(energy,control.drift_weights)/sum(control.drift_weights);
fprintf('Bare and unity transfer: %.9g %.9g\n',bare(1),unity(1));
fprintf('Linked transfer before and after optimisation: %.9g %.9g\n',before(1),after(1));
fprintf('Delivered angular RF energy: %.9g rad^2/s\n',avg_energy);
fprintf('Directional derivative, analytic and finite difference: %.9g %.9g\n',slope,fd);

end


