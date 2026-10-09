% Designs a toroid excitation pulse under Skinner et al. geometry
% Syntax:
%
%    [result,fig]=toroid_oc_optim()
%
% Parameters:
%
%    No input parameters; the RF and offset ensemble is set below
%
% Outputs:
%
%    result - optimised RF waveform and physical response maps
%
%    fig    - radius-offset response and hard-pulse comparison
%
% Source: Skinner et al., J. Magn. Reson. 209, 282-290 (2011).
% DOI: 10.1016/j.jmr.2011.01.026, Eqs. (1),(3),(4).
% This is a new Spinach optimisation under the published probe
% geometry, not the article's original optimised waveform.
%
function [result,fig]=toroid_oc_optim()

% Magnetic field and isotope
sys.magnet=4.7;
sys.isotopes={'1H'};

% Chemical shift, ppm
inter.zeeman.scalar={0};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation={'none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Normalised initial and target states
rho_z=state(spin_system,'Lz',1);
rho_z=rho_z/norm(full(rho_z),2);
rho_x=state(spin_system,'Lx',1);
rho_x=rho_x/norm(full(rho_x),2);

% Control and offset operators
lx=operator(spin_system,'Lx',1);
ly=operator(spin_system,'Ly',1);
lz=operator(spin_system,'Lz',1);

% Drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Sample physical radius uniformly, not the induced RF frequency
radii_m=linspace(1e-3,6e-3,11);

% Control parameters
control.isotopes={'1H'};
control.channels=[1;1];
control.drifts={{H}};
control.operators={lx,ly};
control.off_ops={lz};
control.offsets={linspace(-1.5e3,1.5e3,9)};
control.pwr_levels=2*pi*25e3*6e-3./radii_m;
control.rho_init={rho_z};
control.rho_targ={rho_x};
control.pulse_dt=0.5e-6*ones(1,44);
control.penalties={'SNSA'};
control.p_weights=100;
control.method='lbfgs';
control.max_iter=160;
control.plotting={};

% Reproducible random initial guess
rng(1);
guess=randn(2,numel(control.pulse_dt))/10;

% Spinach optimal-control housekeeping
spin_system=optimcon(spin_system,control);

% Run the optimisation
xy_profile=fmaxnewton(spin_system,@grape_xy,guess);

% Enforce the 25 kHz outer-wall cap
peak_ratio=max(hypot(xy_profile(1,:),xy_profile(2,:)));
rf_hz=25e3*xy_profile/max(1,peak_ratio);

% Verify the shaped pulse and hard benchmark on the same dense grid
check_offsets=linspace(-1.5e3,1.5e3,51);
check_radii=linspace(1e-3,6e-3,51);

% Verify the shaped pulse
[result.response,fig]=toroid_response_map(rf_hz,control.pulse_dt,...
                                          check_offsets,check_radii);

% Simulate the paper's hard-pulse benchmark
[result.baseline,baseline_fig]=toroid_response_map([0;25e3],3.74e-6,...
                                                    check_offsets,check_radii);
close(baseline_fig);

% Package the waveform and optimisation ensemble
result.rf_hz=rf_hz;
result.pulse_dt=control.pulse_dt;
result.radius_m=radii_m;
result.offset_hz=control.offsets{1};
result.peak_hz=max(hypot(rf_hz(1,:),rf_hz(2,:)));

% Compare the new Spinach design with the paper's hard-pulse benchmark
figure(fig);
nexttile(2); hold on;
plot(result.baseline.offset_hz/1e3,result.baseline.detected,'--');
klegend({'Spinach OC','hard benchmark'},'Location','southwest');
ylim([0 1]);

end


