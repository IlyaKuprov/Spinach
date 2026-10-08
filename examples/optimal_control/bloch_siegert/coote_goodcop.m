% Reproduction of the GOODCOP pulse design logic from
% (https://doi.org/10.1038/s41467-018-05400-4) with
% Bloch-Siegert corrections enabled in the optimiser
% and simulator. The pulse enforces contracted-time
% C-alpha evolution while inverting CO. The duration,
% the RF ceiling, the time contraction factor, and the
% two ensemble bands are those of Table 1 and the text
% of the paper; the 53.2 ppm carrier is the one stated
% in Supplementary Figure 5.
%
% Calculation time: minutes
%
% aditya.dev@weizmann.ac.il

function coote_goodcop()

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
alpha_scale=0.90;
pulse_dur=150e-6;
max_rf_hz=15e3;

% Offset grids from the paper
ca_ppm=linspace(35,75,100);
co_ppm=linspace(165,185,30);

% Convert to offset frequencies
ca_hz=ppm2hz(ca_ppm-carrier_ppm,sys.magnet,'13C');
co_hz=ppm2hz(co_ppm-carrier_ppm,sys.magnet,'13C');
all_hz=[ca_hz co_hz];

% Build state pairs correlated to offset ensemble
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

% Set up optimal control problem
control.isotopes={'13C'};           
control.channels=[1; 1];              
control.drifts={{D}};
control.operators={Lx,Ly};
control.rho_init=rho_init;
control.rho_targ=rho_targ;
control.pwr_levels=2*pi*max_rf_hz;
control.pulse_dt=(pulse_dur/75)*ones(1,75);
control.offsets={all_hz};
control.off_ops={Lz};
control.method='lbfgs';
control.max_iter=200;
control.ens_corrs={'rho_ens'};

% Enable BSS corrections
control.bsiegert=true();

% Spinach housekeeping 
spin_system_bs=optimcon(spin_system,control);

% Optimise the waveform
guess=randn(2,numel(control.pulse_dt))/8;
pulse_bs=fmaxnewton(spin_system_bs,@grape_xy,guess);
pulse_bs=control.pwr_levels*pulse_bs;
pulse_bs=mat2cell(pulse_bs,[1 1]);

% Disable BSS corrections
control.bsiegert=false();

% Spinach housekeeping 
spin_system_nobs=optimcon(spin_system,control);

% Optimise the waveform
pulse_nobs=fmaxnewton(spin_system_nobs,@grape_xy,guess);
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
    [controls_aug,pulse_aug]=bloch_siegert(spin_system_bs,{Lx,Ly},pulse_bs);

    % Run with the Bloch-Siegert physics present
    rho=shaped_pulse_xy(spin_system_bs,Hn,controls_aug,pulse_aug,...
                        control.pulse_dt,Iz,'expv-pwc');

    % Get the fidelity
    mz_bs(n)=real(Iz'*rho);

    % Add Bloch-siegert virtual control channels to unadapted pulse
    [controls_aug,pulse_aug]=bloch_siegert(spin_system_bs,{Lx,Ly},pulse_nobs);

    % Run with the Bloch-Siegert physics present
    rho=shaped_pulse_xy(spin_system_bs,Hn,controls_aug,pulse_aug,...
                        control.pulse_dt,Iz,'expv-pwc');

    % Get the fidelity
    mz_nobs(n)=real(Iz'*rho);

end

% Plot the inversion profile
kfigure;
plot(eval_ppm,mz_bs,'LineWidth',1.5); hold on;
plot(eval_ppm,mz_nobs,'LineWidth',1.5); kgrid;
kxlabel('$^{13}$C chemical shift / ppm');
kylabel('Final M$_Z$'); ktitle('GOODCOP profile ');
legend({'BSS on','BSS off'},'Location','Best');
set(gca,'XDir','Reverse');

end


