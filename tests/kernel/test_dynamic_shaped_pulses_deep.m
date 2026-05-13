% Tests dynamic shaped-pulse propagation paths. Syntax:
%
%               result=test_dynamic_shaped_pulses_deep()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks constant-generator reductions of shaped_pulse_xy
% product-quadrature/method combinations and shaped_pulse_af Fokker-Planck
% propagation paths.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_shaped_pulses_deep()

% Announce the test target
fprintf('TESTING: Dynamic shaped pulse propagation paths\n');

% State the dynamic shaped-pulse target of the test
result=new_test_result('kernel/dynamic_shaped_pulses_deep',...
                       'Dynamic shaped pulse propagation paths',...
                       'shaped-pulse helpers must reduce to exact constant-generator propagation when their controls are constant.');

% Build a one-proton Liouville-space spin system
spin_system=local_liouv_system(0);

% Get compact Liouville-space controls
Lx=operator(spin_system,'Lx',1);
Ly=operator(spin_system,'Ly',1);
Lz=operator(spin_system,'Lz',1);
rho=state(spin_system,'Lz',1);

% Define a constant Cartesian RF generator over two slices
slice_durs=[8e-5 13e-5];
amp_x=2*pi*430;
amp_y=-2*pi*170;
drift=2*pi*35*Lz;
generator=drift+amp_x*Lx+amp_y*Ly;
ref_rho=step(spin_system,generator,rho,sum(slice_durs));
ref_prop=propagator(spin_system,generator,sum(slice_durs));

% Exercise the Krylov piecewise-constant path
[rho_obs,traj_obs]=shaped_pulse_xy(spin_system,drift,{Lx,Ly},...
                                   {[amp_x amp_x],[amp_y amp_y]},...
                                   slice_durs,rho,'expv-pwc');
result=test_close(result,'shaped_pulse_xy expv-pwc final',rho_obs,ref_rho,1e-10,1e-10,...
                  'constant PWC controls reduce to one exact Krylov step over the total duration');
result=test_close(result,'shaped_pulse_xy expv-pwc initial trajectory',traj_obs{1},rho,1e-15,1e-15,...
                  'the first trajectory point is the supplied initial state');
result=test_close(result,'shaped_pulse_xy expv-pwc final trajectory',traj_obs{end},ref_rho,1e-10,1e-10,...
                  'the final trajectory point is the propagated final state');

% Exercise the Krylov piecewise-linear path
[rho_obs,traj_obs]=shaped_pulse_xy(spin_system,drift,{Lx,Ly},...
                                   {amp_x*ones(1,3),amp_y*ones(1,3)},...
                                   slice_durs,rho,'expv-pwl');
result=test_close(result,'shaped_pulse_xy expv-pwl final',rho_obs,ref_rho,1e-10,1e-10,...
                  'constant PWL edges reduce to the same exact generator');
result=test_close(result,'shaped_pulse_xy expv-pwl final trajectory',traj_obs{end},ref_rho,1e-10,1e-10,...
                  'PWL trajectory storage returns the final propagated state');

% Exercise explicit exponential product quadratures and propagator output
methods={'expm-pwc','expm-pwl','evol-pwc','evol-pwl'};
for n=1:numel(methods)

    % Select the matching amplitude table shape
    if strcmp(methods{n}(6:8),'pwc')
        amplitudes={[amp_x amp_x],[amp_y amp_y]};
    else
        amplitudes={amp_x*ones(1,3),amp_y*ones(1,3)};
    end

    % Run the shaped pulse with propagator output
    [rho_obs,traj_obs,prop_obs]=shaped_pulse_xy(spin_system,drift,{Lx,Ly},...
                                                amplitudes,slice_durs,rho,methods{n});

    % Check final state and propagator consistency
    result=test_close(result,['shaped_pulse_xy ' methods{n} ' final'],rho_obs,ref_rho,1e-10,1e-10,...
                      'matrix-exponential and evolution paths reduce to exact constant-generator propagation');
    result=test_close(result,['shaped_pulse_xy ' methods{n} ' propagator'],prop_obs,ref_prop,1e-10,1e-10,...
                      'propagator output is the ordered product of constant-slice propagators');
    result=test_close(result,['shaped_pulse_xy ' methods{n} ' trajectory'],traj_obs{end},ref_rho,1e-10,1e-10,...
                      'trajectory storage returns the final propagated state');
end

% Define a constant amplitude-frequency pulse
rf_phi=pi/7;
rf_amp=sqrt(amp_x^2+amp_y^2);
af_generator=rf_amp*(cos(rf_phi)*Lx+sin(rf_phi)*Ly);
af_durs=[6e-5 7e-5 5e-5];
ref_rho=step(spin_system,af_generator,rho,sum(af_durs));
ref_prop=propagator(spin_system,af_generator,sum(af_durs));

% Exercise all shaped_pulse_af propagation methods
methods={'expv','expm','evolution'};
for n=1:numel(methods)

    % Run the Fokker-Planck shaped pulse
    if strcmp(methods{n},'expm')
        [rho_obs,traj_obs,prop_obs]=shaped_pulse_af(spin_system,0*Lx,Lx,Ly,rho,...
                                                    zeros(size(af_durs)),...
                                                    rf_amp*ones(size(af_durs)),...
                                                    af_durs,rf_phi,1,methods{n});
        result=test_close(result,'shaped_pulse_af expm propagator',prop_obs,ref_prop,1e-10,1e-10,...
                          'zero frequency offset and constant RF amplitude project to the exact Liouville propagator');
    else
        [rho_obs,traj_obs]=shaped_pulse_af(spin_system,0*Lx,Lx,Ly,rho,...
                                           zeros(size(af_durs)),...
                                           rf_amp*ones(size(af_durs)),...
                                           af_durs,rf_phi,1,methods{n});
    end

    % Check the folded final state and trajectory
    result=test_close(result,['shaped_pulse_af ' methods{n} ' final'],rho_obs,ref_rho,1e-10,1e-10,...
                      'constant amplitude-frequency controls fold back to exact Cartesian propagation');
    result=test_close(result,['shaped_pulse_af ' methods{n} ' trajectory'],traj_obs{end},ref_rho,1e-10,1e-10,...
                      'the folded trajectory final point matches the returned final state');
end

end

% Local quiet one-spin Liouville-space test system
function spin_system=local_liouv_system(magnet)

% Specify the spin system
sys.magnet=magnet;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};

% Specify a full Liouville basis
bas.formalism='zeeman-liouv';
bas.approximation='none';

% Build the quiet regression-test system
spin_system=test_spin_system(sys,inter,bas);

end


