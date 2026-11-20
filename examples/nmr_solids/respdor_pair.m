% Compares R- and S-RESPDOR dephasing in a model 13C-14N pair
% at matched field, MAS rate, offsets, and saturation RF amplitude.
% Uses the Grenoble R and S pulse schedules with common interactions
% and echo RF. A reduced quadrupolar coupling and a small powder
% grid keep this first-order high-field demonstration short.
%
% talos@spindynamics.org

function respdor_pair()

% Model carbon-nitrogen pair
sys.magnet=2*pi*600e6/spin('1H');
sys.isotopes={'13C','14N'};
inter.coupling.eigs={[] [-1000 -1000 2000]; [] []};
inter.coupling.euler={[] deg2rad([0 30 0]); [] []};
inter.coupling.matrix=cell(2,2);
inter.coupling.matrix{2,2}=eeqq2nqi(100e3,0,1,[0 0 0]);

% Basis and local workers
bas.formalism='sphten-liouv';
bas.approximation={'none'};
sys.parallel={'processes',4};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Shared experimental settings
parameters.spins=sys.isotopes;
parameters.rate=10000;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.offset=[0 0];
parameters.max_rank=12;
parameters.grid='leb_2ang_rank_11';
parameters.echo_rf=100e3;
parameters.sat_rf=40e3;

% R-RESPDOR with longitudinal preparation and detection
parameters.sequence='R';
parameters.ncycles=1:12;
parameters.sat_time=2/parameters.rate;
parameters.rho0=state(spin_system,'Lz','13C');
parameters.coil=-coil_state(spin_system,'Lz','13C');
curve_r=singlerot(spin_system,@respdor,parameters,'nmr');
time_r=4*parameters.ncycles/parameters.rate;

% S-RESPDOR with transverse preparation and detection
parameters.sequence='S';
parameters.ncycles=0:8;
parameters.sat_time=1/parameters.rate;
parameters.rho0=state(spin_system,'Lx','13C');
parameters.coil=coil_state(spin_system,'L+','13C');
curve_s=singlerot(spin_system,@respdor,parameters,'nmr');
time_s=6*parameters.ncycles/parameters.rate;

% Normalise only after powder averaging and omit unresolved references
contrast_r=nan(size(time_r));
contrast_s=nan(size(time_s));
valid_r=abs(curve_r(1,:))>spin_system.tols.liouv_zero;
valid_s=abs(curve_s(1,:))>spin_system.tols.liouv_zero;
contrast_r(valid_r)=real(curve_r(3,valid_r)./curve_r(1,valid_r));
contrast_s(valid_s)=real(curve_s(3,valid_s)./curve_s(1,valid_s));

% Plot contrast against the combined duration of the two recoupling arms
kfigure(); plot(1e3*time_r,contrast_r,'o-',1e3*time_s,contrast_s,'s-');
kxlabel('Total recoupling time, ms');
kylabel('$\Delta S/S_0$'); klegend({'R-RESPDOR','S-RESPDOR'}); kgrid;

end


