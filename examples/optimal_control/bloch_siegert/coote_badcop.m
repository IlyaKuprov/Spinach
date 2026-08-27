% Reproduction of BADCOP-style selective decoupling logic from Coote et al.
% with Bloch-Siegert corrections enabled in the optimiser and simulator
% BADCOP1, BADCOP2, and BADCOP3 are designed and validated
%
% aditya.dev@weizmann.ac.il

function coote_badcop()

% Magnetic field corresponding to 800 MHz 1H
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

% Relevant operators and states
Lx=operator(spin_system,'Lx','13C');
Ly=operator(spin_system,'Ly','13C');
Lz=operator(spin_system,'Lz','13C');
Ix=state(spin_system,'Lx','13C'); Ix=Ix/norm(full(Ix),2);
Iy=state(spin_system,'Ly','13C'); Iy=Iy/norm(full(Iy),2);
Iz=state(spin_system,'Lz','13C'); Iz=Iz/norm(full(Iz),2);

% Drift Hamiltonian
D=hamiltonian(assume(spin_system,'nmr'));

% Shared paper parameters
carrier_ppm=100;
alpha_scale=0.91;
pulse_dur=1e-3;
ca_ppm=linspace(40,72,100);
co_ppm=linspace(165,185,30);
ca_hz=ppm2hz(sys.magnet,'13C',ca_ppm,carrier_ppm);
co_hz=ppm2hz(sys.magnet,'13C',co_ppm,carrier_ppm);

% Build all three variants from Table 1
variants={...
    struct('name','BADCOP1','rf_hz',5.94e3,'cb_inv_ppm',[5 37]),...
    struct('name','BADCOP2','rf_hz',4.87e3,'cb_inv_ppm',[28 35]),...
    struct('name','BADCOP3','rf_hz',7.22e3,'cb_inv_ppm',[10 45])};

% Design and evaluate each variant
for k=1:numel(variants)

    % Build C-beta inversion and preservation grids
    cb_inv_ppm=linspace(variants{k}.cb_inv_ppm(1),variants{k}.cb_inv_ppm(2),60);
    cb_prs_ppm=linspace(5,80,80);
    mask=(cb_prs_ppm<variants{k}.cb_inv_ppm(1))|(cb_prs_ppm>variants{k}.cb_inv_ppm(2));
    cb_prs_ppm=cb_prs_ppm(mask);
    cb_inv_hz=ppm2hz(sys.magnet,'13C',cb_inv_ppm,carrier_ppm);
    cb_prs_hz=ppm2hz(sys.magnet,'13C',cb_prs_ppm,carrier_ppm);

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
    control.max_iter=120;
    control.ens_corrs={'rho_ens'};
    control.plotting={};

    % Enable Bloch-Siegert corrections
    control.bsiegert=true();

    % Spinach housekeeping for Bloch-Siegert case
    ss_bs=optimcon(spin_system,control);

    % Optimise the Bloch-Siegert waveform
    guess=randn(2,numel(control.pulse_dt))/8;
    pulse_bs=fmaxnewton(ss_bs,@grape_xy,guess);
    pulse_bs=control.pwr_levels*pulse_bs;
    pulse_bs=mat2cell(pulse_bs,[1 1]);

    % Disable Bloch-Siegert corrections
    control.bsiegert=false();

    % Spinach housekeeping for no Bloch-Siegert case
    ss_nobs=optimcon(spin_system,control);

    % Optimise the no Bloch-Siegert waveform
    pulse_nobs=fmaxnewton(ss_nobs,@grape_xy,guess);
    pulse_nobs=control.pwr_levels*pulse_nobs;
    pulse_nobs=mat2cell(pulse_nobs,[1 1]);

    % Validate profile on a dense ppm grid
    eval_ppm=linspace(0,200,251);
    eval_hz=ppm2hz(sys.magnet,'13C',eval_ppm,carrier_ppm);
    mz_bs=zeros(size(eval_hz));
    mz_nobs=zeros(size(eval_hz));
    for n=1:numel(eval_hz)
        Hn=D+2*pi*eval_hz(n)*Lz;


        [controls_aug,pulse_aug]=bloch_siegert(ss_bs,{Lx,Ly},...
                                                     pulse_bs);
        rho=shaped_pulse_xy(ss_bs,Hn,controls_aug,pulse_aug,...
                            control.pulse_dt,Iz,'expv-pwc');
        mz_bs(n)=real(Iz'*rho);
        [controls_aug,pulse_aug]=bloch_siegert(ss_bs,{Lx,Ly},...
                                                     pulse_nobs);
        rho=shaped_pulse_xy(ss_bs,Hn,controls_aug,pulse_aug,...
                            control.pulse_dt,Iz,'expv-pwc');
        mz_nobs(n)=real(Iz'*rho);
    end

    % Plot the inversion profile
    figure;
    plot(eval_ppm,mz_bs,'LineWidth',1.5); hold on;
    plot(eval_ppm,mz_nobs,'LineWidth',1.5);
    xline(variants{k}.cb_inv_ppm(1),'k--','LineWidth',1.0);
    xline(variants{k}.cb_inv_ppm(2),'k--','LineWidth',1.0);
    grid on;
    xlabel('^{13}C chemical shift / ppm');
    ylabel('Final M_z');
    title([variants{k}.name ' profile with and without Bloch-Siegert']);
    legend({'With Bloch-Siegert','No Bloch-Siegert'},'Location','Best');

end

end

function off_hz=ppm2hz(magnet,isotope,ppm_grid,carrier_ppm)

% Convert ppm values into offset frequencies in Hz
frq_hz=abs(spin(isotope)*magnet/(2*pi));
off_hz=(ppm_grid-carrier_ppm)*1e-6*frq_hz;

end


