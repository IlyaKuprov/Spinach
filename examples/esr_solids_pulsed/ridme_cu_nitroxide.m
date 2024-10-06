% RIDME on a Cu(II)-NO two electron system at Q-band.
%
% The numerical calculation is done by brute-force time propaga-
% tion and numerical powder averaging in Liouville space, inclu-
% ding g-factor orientation effects on the dipolar coupling.
%
% Relaxation is included as extended T1/T2 theory.
%
% The analytical calculation is done for isotropic parts of the
% electron g-factors.
%
% Calculation time: seconds.
%
% alice.bowen@chem.ox.ac.uk
% i.kuprov@soton.ac.uk

function ridme_cu_nitroxide()

% Spin system parameters
sys.magnet=1.249;
sys.isotopes={'E','E'};                     % Spin 1 is Cu, Spin 2 is NO
inter.zeeman.eigs={[2.056, 2.056, 2.205];   % Copper(II)
                   [2.009, 2.006, 2.003]};  % Nitroxide
inter.zeeman.euler={[0 0 0]; [0 0 0]};
inter.coordinates={[43 0 0]; [0 0 0]};      % Angstrom

% Relaxation theory
inter.relaxation={'t1_t2'};
inter.r1_rates={1/(35e-6)  1/(2e-3)};     % rates in Hz
inter.r2_rates={1/(1.5e-6) 1/(1.3e-6)};   % rates in Hz
inter.rlx_keep='labframe';
inter.equilibrium='zero';

% Formalism
bas.formalism='sphten-liouv';
bas.approximation='none';

% Disable trajectory level SSR algorithms
sys.disable={'trajlevel'};
               
% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.rho0=state(spin_system,'Lz','E');  % initial condition
parameters.probe_spin=2;                      % probe spin
parameters.stepsize=16e-9;                    % timestep
parameters.nsteps=[25 188];                   % number of steps in tau1 and tau2
parameters.tmix=35e-6;                        % mixing time 
parameters.spins={'E'};
parameters.grid='rep_2ang_800pts_sph';

% Simulation
answer=powder(spin_system,@ridme,parameters,'esr');

% Phase cycle processing
answer.ridme_trace.real=answer.pxpxpx.real+answer.pypypx.real+ ...
                        answer.mxmxpx.real+answer.mymypx.real;
answer.ridme_trace.imag=answer.pxpxpx.imag+answer.pypypx.imag+ ...
                        answer.mxmxpx.imag+answer.mymypx.imag;

% Axis ticks and figure frame
x_axis=linspace((-parameters.stepsize*parameters.nsteps(1)),...
                 (parameters.stepsize*parameters.nsteps(2)),...
                  parameters.nsteps(1)+parameters.nsteps(2)+1);
figure(); scale_figure([1.50 0.75]);

% Real parts
subplot(1,2,1); hold on; kgrid; box on;
plot(1e6*x_axis,real(answer.pxpxpx.real));
plot(1e6*x_axis,real(answer.pypypx.real));
plot(1e6*x_axis,real(answer.mxmxpx.real));
plot(1e6*x_axis,real(answer.mymypx.real));
plot(1e6*x_axis,real(answer.ridme_trace.real));
kxlabel('time, $\mu$s'); ktitle('real parts'); xlim tight;
klegend({'PxPxPx','PyPyPx','MxMxPx','MyMyPx','RIDME'});

% Imaginary parts
subplot(1,2,2); hold on; kgrid; box on;
plot(1e6*x_axis,real(answer.pxpxpx.imag));
plot(1e6*x_axis,real(answer.pypypx.imag));
plot(1e6*x_axis,real(answer.mxmxpx.imag));
plot(1e6*x_axis,real(answer.mymypx.imag));
plot(1e6*x_axis,real(answer.ridme_trace.imag));
kxlabel('time, $\mu$s'); ktitle('imag parts'); xlim tight;
klegend({'PxPxPx','PyPyPx','MxMxPx','MyMyPx','RIDME'});

end

