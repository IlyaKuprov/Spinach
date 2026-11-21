% Finite-gradient CRAZED echo build-up in CuSO4-doped water, Fig. 2 of
% Barros et al., https://doi.org/10.1063/1.3116107. The Fig. 1b sequence
% is repeated for each echo time and each of four excitation phases.
% The receiver is calibrated by a separate 90-degree, gradient-free
% acquisition; there is no fit to the digitised experimental points.
%
% Example: [residuals,peaks]=ddf_barros_water(2,16,0.002,8);
% periods is the number of modulation wavelengths in the reduced 1D
% grating (zero instead selects the full 8 mm 3D sphere). Use a positive
% integer points_per_wave, a positive max_step in seconds, and a positive
% integer workers count. The measured magnetisation is in units of M0.
% The full sphere at resolved pitch requires a large-memory solver.
%
% residuals contains the pointwise Fig. 2 comparisons in M0 units; peaks
% compares sampled experimental and model maxima and their times, with
% signed model-minus-experiment differences. traces contains the four
% raw complex shot signals and the complex 90-degree reference, indexed
% by echo time (rows), gradient and orientation (columns), and shot (pages). The RF widths and phase cycle
% phases are not specified in the paper; short rectangular RF pulses and
% an equivalent ideal second-order coherence selection are used here.

function [residuals,peaks,traces]=ddf_barros_water(periods,points_per_wave,max_step,workers)

% Doped water parameters and experimental settings from Section 6 and Fig. 2
sys.magnet=9.4;
sys.isotopes={'1H'};
sys.parallel={'processes',workers};
sys.output='hush';
inter.zeeman.scalar={0};
t_one=0.288;
t_two=0.250;
tau_d=0.098;
diffusivity=2.15e-9;
diameter=8e-3;
delta=2e-3;
t_first=5e-3;
rf_width=10e-6;
grad_amp=[0.05 0.10 0.15 0.20];
inv_weight=[6.58 1.64 0.73 0.41];
direction=repelem([0 0 1; 1 0 0],4,1);
phase=exp(1i*pi*(0:3)/2);
weights=[1 -1 1 -1]/4;

% Single-proton resonant rotating frame and complete Liouville basis
bas.formalism='sphten-liouv';
bas.approximation={'none'};
spin_system=basis(create(sys,inter),bas);
gamma=spin('1H');
magn_eq=1/(spin_system.tols.mu0*gamma*tau_d);

% Affine Bloch recovery represented linearly by retaining the identity
rho_id=state(spin_system,'E','1H');
rho_z=state(spin_system,'Lz','1H');
detect=coil_state(spin_system,'L+','1H');
probe_z=coil_state(spin_system,'Lz','1H');
proj_id=rho_id'/(rho_id'*rho_id);
proj_z=rho_z'/(rho_z'*rho_z);
dim=numel(rho_id);
relax=-(speye(dim)-rho_id*proj_id-rho_z*proj_z)/t_two+...
      rho_z*(proj_id-proj_z)/t_one;
params.groups={1};
params.frames=1;
params.scale=magn_eq/real(probe_z'*rho_z);

% Select simulated echo times and order the eight experimental series
time_axis=[18:20:598 612]'/1000;
signal=zeros(numel(time_axis),8);
first_order=zeros(size(signal));
shots=complex(zeros([size(signal) 4]));
references=complex(zeros(1,8));
curve_grad=repmat(grad_amp,1,2);
curve_weight=repmat(inv_weight,1,2);
curve_panel=repelem(["a" "b"],4);

% Keep only one large sphere model in memory at a time
curve_workers=workers*(periods>0);
parfor (k=1:8,curve_workers)

    % Set the physical extent and pitch from the gradient modulation wavelength
    axis_dir=direction(k,:);
    wave_num=1/sqrt(tau_d*diffusivity*curve_weight(k));
    wavelength=2*pi/wave_num;
    parameters=params;
    if periods>0

        % A plane grating with a spherical mean field
        parameters.npts=1+(periods*points_per_wave-1)*axis_dir;
        spacing=repmat(wavelength/points_per_wave,1,3);
        parameters.dims=parameters.npts.*spacing;
        mask=true(parameters.npts);
        parameters.field=@(magn)ddf_field(magn,axis_dir,1/3);

    else

        % Voxelised sphere with vacuum outside and reflecting diffusion inside
        parameters.npts=repmat(ceil(diameter*points_per_wave/wavelength),1,3);
        parameters.dims=repmat(diameter,1,3);
        spacing=parameters.dims./parameters.npts;
        axis_pos=((1:parameters.npts(1))-(parameters.npts(1)+1)/2)*spacing(1);
        [xx,yy,zz]=ndgrid(axis_pos,axis_pos,axis_pos);
        mask=(xx.^2+yy.^2+zz.^2)<=(diameter/2)^2;
        symbol=ddf_bounded(parameters.npts,spacing);
        parameters.field=@(magn)ddf_cartesian(magn,parameters.npts,symbol);

    end

    % Prepare the volume-averaged receiver and fixed diffusion and relaxation
    parameters.spc_dim=prod(parameters.npts);
    parameters.spn_dim=dim;
    parameters.rho0=kron(double(mask(:)),rho_id+rho_z);
    parameters.coil=kron(double(mask(:))/nnz(mask),detect);
    diffusion=ddf_diffusion(mask,spacing,diffusivity);
    if periods>0

        % Close the diffusion line periodically across one end face
        seam=sparse([1 parameters.spc_dim],[1 1],[1 -1],parameters.spc_dim,1);
        diffusion=diffusion-(diffusivity/spacing(1)^2)*(seam*seam.');

    end
    parameters.transport=1i*kron(diffusion,speye(dim));
    R=kron(spdiags(double(mask(:)),0,parameters.spc_dim,parameters.spc_dim),relax);
    empty=sparse(parameters.spc_dim*dim,parameters.spc_dim*dim);
    G=cellfun(@inflate,g2fplanck(spin_system,parameters),'UniformOutput',false);

    % Calibrate the fixed receiver on a hard 90-degree readout with no feedback
    p_ref=parameters;
    p_ref.scale=0;
    p_ref.duration=rf_width;
    p_ref.rf=pi/(2*rf_width);
    p_ref.gradient=zeros(1,3);
    p_ref.steps=1;
    p_ref.acquire=true;
    reference=ddf_sequence(spin_system,p_ref,empty,empty,empty,G,empty);
    model={parameters,empty,R,G,wave_num*axis_dir/(gamma*delta),reference};

    % Acquire all four completed transients before forming each echo
    [signal(:,k),shots(:,k,:)]=echo_curve(spin_system,model,time_axis,delta,...
                                         t_first,rf_width,max_step,phase,weights);
    references(k)=reference;

    % Evaluate the published first-order, instantaneous-refocusing approximation
    rate=2*diffusivity*wave_num^2+1/t_one;
    angular=(3*axis_dir(3)^2-1)/2;
    first_order(:,k)=abs(angular)/(4*tau_d)*...
                     exp(-2*(diffusivity*wave_num^2*(t_first-delta/3)+t_first/t_two))*...
                     exp(-time_axis/t_two).*(-expm1(-rate*time_axis))/rate;

end

% Preserve the uncombined complex acquisitions and their calibration
traces.t2_s=time_axis;
traces.panel=curve_panel;
traces.gradient_T_m=curve_grad;
traces.shots=shots;
traces.references=references;
traces.weights=weights;

% Compare at digitised times and report signed peak and pointwise discrepancies
data_file=fullfile(fileparts(mfilename('fullpath')),'barros_fig2.csv');
residuals=readtable(data_file,'TextType','string');
residuals.model=zeros(height(residuals),1);
residuals.eq_one=zeros(height(residuals),1);
summary=zeros(8,12);
for k=1:8

    % Interpolate the two curves to the times of the open-circle data
    rows=(residuals.panel==curve_panel(k))&(residuals.gradient_T_per_m==curve_grad(k));
    times=residuals.t2_ms(rows);
    model=interp1(time_axis,signal(:,k),times/1000,'pchip','extrap');
    eq_one=interp1(time_axis,first_order(:,k),times/1000,'pchip','extrap');
    residuals.model(rows)=model;
    residuals.eq_one(rows)=eq_one;
    observed=residuals.signal(rows);
    [exp_peak,j]=max(observed);
    amplitude_unc=residuals.signal_unc(rows);
    time_unc=residuals.t2_unc_ms(rows);
    low=max(observed-amplitude_unc);
    high=max(observed+amplitude_unc);
    possible=(observed+amplitude_unc)>=low;
    [sim_peak,n]=max(model);
    summary(k,:)=[sqrt(mean((model-observed).^2)),...
                  sqrt(mean((eq_one-observed).^2)),...
                  exp_peak sim_peak sim_peak-exp_peak,...
                  times(j) times(n) times(n)-times(j) low high,...
                  min(times(possible)-time_unc(possible)),...
                  max(times(possible)+time_unc(possible))];

end

% Store signed residuals and display the eight separate comparisons
residuals.error=residuals.model-residuals.signal;
residuals.eq_one_error=residuals.eq_one-residuals.signal;
peaks=[table(curve_panel.',curve_grad.','VariableNames',{'panel','G_T_per_m'}),...
       array2table(summary,'VariableNames',{'rmse','eq_one_rmse','exp_peak','sim_peak',...
                                           'peak_error','exp_ms','sim_ms','peak_ms_error',...
                                           'exp_peak_low','exp_peak_high',...
                                           'exp_ms_low','exp_ms_high'})];
disp(peaks);

% Plot experimental points, simulations, and published Eq. 1 on identical axes
kfigure();
colours=lines(4);
tiledlayout(1,2);
for j=1:2

    % Each panel holds one gradient orientation relative to B0
    nexttile();
    hold on
    for k=1:4

        % Solid: production solver; dashed: Eq. 1; circles: experiment
        n=(j-1)*4+k;
        rows=(residuals.panel==curve_panel(n))&(residuals.gradient_T_per_m==curve_grad(n));
        plot(time_axis,signal(:,n),'Color',colours(k,:),'DisplayName',num2str(curve_grad(n)));
        plot(time_axis,first_order(:,n),'--','Color',colours(k,:),'HandleVisibility','off');
        plot(residuals.t2_ms(rows)/1000,residuals.signal(rows),'o','Color',colours(k,:),...
             'HandleVisibility','off');

    end
    kxlabel('$t_2$, s');
    kylabel('echo amplitude / $M_0$');
    orientations={'parallel to $B_0$','perpendicular to $B_0$'};
    ktitle(orientations{j});
    leg=klegend('Location','northeast');
    leg.Title.String='G, T/m';
    ylim([0 0.14]);
    kgrid

end

end

% Simulate and cycle each echo independently
function [signal,shots]=echo_curve(spin_system,model,time_axis,delta,t_first,...
                                              rf_width,max_step,phase,weights)

% Initialise the echo curve
signal=zeros(size(time_axis));
shots=complex(zeros(numel(time_axis),1,4));
for n=1:numel(time_axis)

    % Set the refocusing pulse at the midpoint of the two equal tau intervals
    t_echo=time_axis(n);
    tau=(t_echo-2*delta-2*t_first)/2;
    parameters=model{1};
    empty=model{2};
    R=model{3};
    G=model{4};
    gradient=model{5};
    reference=model{6};
    parameters.duration=[rf_width; delta; t_first-delta-rf_width; rf_width;...
                         2*delta; 2*t_first+tau-rf_width/2; rf_width; tau-rf_width/2];
    parameters.gradient=zeros(numel(parameters.duration),3);
    parameters.gradient([2 5],:)=[gradient; gradient];
    parameters.steps=max(1,ceil(parameters.duration/max_step));
    parameters.acquire=false(size(parameters.duration));
    parameters.acquire(end)=true;

    % Project onto second-order coherence before taking the magnitude
    fid=zeros(1,4);
    for q=1:4

        % Excitation phase cycle, fixed minus-Y mixing, and fixed X refocusing
        parameters.rf=zeros(size(parameters.duration));
        parameters.rf([1 4 7])=[phase(q)*pi/2; -1i*pi/2; pi]/rf_width;
        fid(q)=ddf_sequence(spin_system,parameters,empty,R,empty,G,empty);

    end
    shots(n,1,:)=fid;
    signal(n)=abs(sum(weights.*fid)/reference);

end

end


