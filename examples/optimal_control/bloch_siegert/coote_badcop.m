% Reproduction of BADCOP-style selective decoupling logic
% from (https://doi.org/10.1038/s41467-018-05400-4) with
% Bloch-Siegert corrections enabled in the optimiser and
% simulator. BADCOP1, BADCOP2, and BADCOP3 are designed
% and validated. Durations, RF ceilings, the contraction
% factor, and the inversion bands are those of Table 1
% and the text of the paper, and the 53.2 ppm carrier is
% the one stated in Supplementary Figure 5. Only BADCOP2
% and BADCOP3 are asked to preserve C-beta magnetisation
% outside their inversion band.
%
% Calculation time: minutes
%
% aditya.dev@weizmann.ac.il

function coote_badcop()

% Magnetic field
sys.magnet=18.8;

% Single-spin carbon model
sys.isotopes={'13C'};
inter.zeeman.scalar={0.0};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get relevant operators
Lx=operator(spin_system,'Lx','13C');
Ly=operator(spin_system,'Ly','13C');
Lz=operator(spin_system,'Lz','13C');

% Get and normalise relevant states
Ix=state(spin_system,'Lx','13C'); Ix=Ix/norm(Ix,2);
Iy=state(spin_system,'Ly','13C'); Iy=Iy/norm(Iy,2);
Iz=state(spin_system,'Lz','13C'); Iz=Iz/norm(Iz,2);

% Drift Hamiltonian
D=hamiltonian(assume(spin_system,'nmr'));

% Paper parameters
carrier_ppm=53.2;
alpha_scale=0.91;
pulse_dur=1e-3;
ca_ppm=linspace(40,72,100);
co_ppm=linspace(165,185,30);
ca_hz=ppm2hz(ca_ppm-carrier_ppm,sys.magnet,'13C');
co_hz=ppm2hz(co_ppm-carrier_ppm,sys.magnet,'13C');

% Build all three variants from Table 1
variants={...
    struct('name','BADCOP1','rf_hz',5.94e3,'cb_inv_ppm',[5 37],'cb_prs',false),...
    struct('name','BADCOP2','rf_hz',4.87e3,'cb_inv_ppm',[28 35],'cb_prs',true),...
    struct('name','BADCOP3','rf_hz',7.22e3,'cb_inv_ppm',[10 45],'cb_prs',true)};

% Design and evaluate each variant
for k=1:numel(variants)

    % Build the C-beta inversion grid
    cb_inv_ppm=linspace(variants{k}.cb_inv_ppm(1),...
                        variants{k}.cb_inv_ppm(2),60);
    cb_inv_hz=ppm2hz(cb_inv_ppm-carrier_ppm,sys.magnet,'13C');

    % Only BADCOP2 and BADCOP3 preserve C-beta outside that band
    if variants{k}.cb_prs
        cb_prs_ppm=linspace(5,80,80);
        mask=(cb_prs_ppm<variants{k}.cb_inv_ppm(1))|...
             (cb_prs_ppm>variants{k}.cb_inv_ppm(2));
        cb_prs_hz=ppm2hz(cb_prs_ppm(mask)-carrier_ppm,sys.magnet,'13C');
    else
        cb_prs_hz=[];
    end

    % Assemble full offset list
    all_hz=[ca_hz co_hz cb_inv_hz cb_prs_hz];

    % Build correlated state-to-state targets
    rho_init=cell(1,numel(all_hz));
    rho_targ=cell(1,numel(all_hz));
    for n=1:numel(ca_hz)
        if mod(n,2)==1
            rho_init{n}=Ix;
        else
            rho_init{n}=Iy;
        end
        rho_targ{n}=step(spin_system,2*pi*ca_hz(n)*Lz,rho_init{n},alpha_scale*pulse_dur);
    end
    for n=1:numel(co_hz)
        idx=numel(ca_hz)+n;
        rho_init{idx}=Iz;
        rho_targ{idx}=-Iz;
    end
    for n=1:numel(cb_inv_hz)
        idx=numel(ca_hz)+numel(co_hz)+n;
        rho_init{idx}=Iz;
        rho_targ{idx}=-Iz;
    end
    for n=1:numel(cb_prs_hz)
        idx=numel(ca_hz)+numel(co_hz)+numel(cb_inv_hz)+n;
        rho_init{idx}=Iz;
        rho_targ{idx}=Iz;
    end

    % Set up optimal control problem
    control.isotopes={'13C'};
    control.channels=[1; 1];
    control.drifts={{D}};
    control.operators={Lx,Ly};
    control.rho_init=rho_init;
    control.rho_targ=rho_targ;
    control.pwr_levels=2*pi*variants{k}.rf_hz;
    control.pulse_dt=(pulse_dur/200)*ones(1,200);
    control.offsets={all_hz};
    control.off_ops={Lz};
    control.method='lbfgs';
    control.max_iter=200;
    control.ens_corrs={'rho_ens'};
    control.plotting={};

    % Enable Bloch-Siegert corrections
    control.bsiegert=true();

    % Spinach housekeeping 
    ss_bs=optimcon(spin_system,control);

    % Optimise the Bloch-Siegert waveform
    guess=randn(2,numel(control.pulse_dt))/8;
    pulse_bs=fmaxnewton(ss_bs,@grape_xy,guess);
    pulse_bs=control.pwr_levels*pulse_bs;
    pulse_bs=mat2cell(pulse_bs,[1 1]);

    % Disable Bloch-Siegert corrections
    control.bsiegert=false();

    % Spinach housekeeping 
    ss_nobs=optimcon(spin_system,control);

    % Optimise the no Bloch-Siegert waveform
    pulse_nobs=fmaxnewton(ss_nobs,@grape_xy,guess);
    pulse_nobs=control.pwr_levels*pulse_nobs;
    pulse_nobs=mat2cell(pulse_nobs,[1 1]);

    % Validate on a dense ppm grid
    eval_ppm=linspace(0,200,251);
    eval_hz=ppm2hz(eval_ppm-carrier_ppm,sys.magnet,'13C');
    mz_bs=zeros(size(eval_hz));
    mz_nobs=zeros(size(eval_hz));
    for n=1:numel(eval_hz)

        % Add the offset term
        Hn=D+2*pi*eval_hz(n)*Lz;

        % Add Bloch-siegert virtual control channels to adapted pulse
        [controls_aug,pulse_aug]=bloch_siegert(ss_bs,{Lx,Ly},pulse_bs);
        
        % Run with the Bloch-Siegert physics present
        rho=shaped_pulse_xy(ss_bs,Hn,controls_aug,pulse_aug,...
                            control.pulse_dt,Iz,'expv-pwc');

        % Get the fidelity
        mz_bs(n)=real(Iz'*rho);

        % Add Bloch-siegert virtual control channels to unadapted pulse
        [controls_aug,pulse_aug]=bloch_siegert(ss_bs,{Lx,Ly},pulse_nobs);

        % Run with the Bloch-Siegert physics present
        rho=shaped_pulse_xy(ss_bs,Hn,controls_aug,pulse_aug,...
                            control.pulse_dt,Iz,'expv-pwc');

        % Get the fidelity
        mz_nobs(n)=real(Iz'*rho);

    end

    % Plot the inversion profile
    kfigure;
    plot(eval_ppm,mz_bs,'LineWidth',1.5); hold on;
    plot(eval_ppm,mz_nobs,'LineWidth',1.5); kgrid
    xline(variants{k}.cb_inv_ppm(1),'k--','LineWidth',1.0);
    xline(variants{k}.cb_inv_ppm(2),'k--','LineWidth',1.0);
    kxlabel('$^{13}$C chemical shift / ppm'); 
    kylabel('Final M$_Z$'); ktitle([variants{k}.name ' profile']);
    klegend({'BSS on','BSS off'},'Location','Best'); 
    set(gca,'XDir','Reverse'); drawnow;

end

end


