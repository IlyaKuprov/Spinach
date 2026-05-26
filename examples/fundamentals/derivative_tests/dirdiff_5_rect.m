% GRAPE Hessian internal consistency test: Newton 
% against Goodwin algorithm.
%
% ilya.kuprov@weizmann.ac.il
% david.goodwin@inano.au.dk

function dirdiff_5_rect()

% Formalisms to test
formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'};

% Loop over formalisms
for n=1:numel(formalisms)

    % Build the derivative-test system
    [spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalisms{n});

    % Define control parameters
    control.drifts={{H}};                           % Drift
    control.operators={Lx,Ly};                      % Controls
    control.rho_init={ Sx Sy Sz};                   % Starting states
    control.rho_targ={-Sz Sy Sx};                   % Target states
    control.pwr_levels=2*pi*linspace(50e3,70e3,10); % Power levels
    control.max_iter=1000;                          % Termination condition
    control.plotting={};                            % Plotting options
    control.integrator='rectangle';                 % Integrator
    control.pulse_dt=12.8e-6*ones(1,5);             % Time interval grid

    % Pick initial guess, phase-modulated GRAPE
    control.amplitudes=ones(1,5); guess=randn(1,5)/3;

    % Get Newton Hessian
    control.method='newton';
    spin_system=optimcon(spin_system,control);
    [~,~,~,newton_hess_ph]=grape_phase(guess,spin_system);

    % Get Goodwin Hessian
    control.method='goodwin';
    spin_system=optimcon(spin_system,control);
    [~,~,~,goodwin_hess_ph]=grape_phase(guess,spin_system);

    % Pick initial guess, XY-modulated GRAPE
    guess=randn(2,5)/3;

    % Get Newton Hessian
    control.method='newton';
    spin_system=optimcon(spin_system,control);
    [~,~,~,newton_hess_xy]=grape_xy(guess,spin_system);

    % Get Goodwin Hessian
    control.method='goodwin';
    spin_system=optimcon(spin_system,control);
    [~,~,~,goodwin_hess_xy]=grape_xy(guess,spin_system);

    % Run the comparisons
    if norm(newton_hess_ph(:)-goodwin_hess_ph(:),1)>1e-6*norm(newton_hess_ph(:),1)
        error([formalisms{n} ' phase Hessian internal consistency test failed.']);
    else
        disp([formalisms{n} ' phase Hessian internal consistency test passed.']);
    end
    if norm(newton_hess_xy(:)-goodwin_hess_xy(:),1)>1e-6*norm(newton_hess_xy(:),1)
        error([formalisms{n} ' Cartesian Hessian internal consistency test failed.']);
    else
        disp([formalisms{n} ' Cartesian internal consistency test passed.']);
    end

end

end
