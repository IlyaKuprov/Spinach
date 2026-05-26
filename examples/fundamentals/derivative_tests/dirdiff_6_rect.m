% GRAPE Hessian test against finite-differenced gradients.
%
% ilya.kuprov@weizmann.ac.il
% david.goodwin@inano.au.dk

function dirdiff_6_rect()

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
    control.method='newton';                        % Optimisation method

    % Spinach housekeeping
    spin_system=optimcon(spin_system,control);

    % Random guess and finite diff increment
    guess=randn(2,5)/3; h=sqrt(eps('double'));

    % Call GRAPE and request analytical Hessian
    [~,~,~,hess_anl]=grape_xy(guess,spin_system);
    hess_anl=squeeze(hess_anl(:,:,1));

    % Leftmost Hessian column
    wave_forw=guess; wave_forw(1)=wave_forw(1)+h;
    wave_back=guess; wave_back(1)=wave_back(1)-h;
    [~,~,grad_forw]=grape_xy(wave_forw,spin_system);
    [~,~,grad_back]=grape_xy(wave_back,spin_system);
    grad_forw=squeeze(grad_forw(:,:,1)); grad_forw=grad_forw(:);
    grad_back=squeeze(grad_back(:,:,1)); grad_back=grad_back(:);
    hess_num=(grad_forw-grad_back)/(2*h);
    if norm(hess_anl(:,1)-hess_num,1)<1e-6*norm(hess_num,1)
        disp([formalisms{n} ' leftmost column test passed']);
    else
        error([formalisms{n} ' leftmost column test failed']);
    end

    % Rightmost Hessian column
    wave_forw=guess; wave_forw(end)=wave_forw(end)+h;
    wave_back=guess; wave_back(end)=wave_back(end)-h;
    [~,~,grad_forw]=grape_xy(wave_forw,spin_system);
    [~,~,grad_back]=grape_xy(wave_back,spin_system);
    grad_forw=squeeze(grad_forw(:,:,1)); grad_forw=grad_forw(:);
    grad_back=squeeze(grad_back(:,:,1)); grad_back=grad_back(:);
    hess_num=(grad_forw-grad_back)/(2*h);
    if norm(hess_anl(:,end)-hess_num,1)<1e-6*norm(hess_num,1)
        disp([formalisms{n} ' rightmost column test passed']);
    else
        error([formalisms{n} ' rightmost column test failed']);
    end

    % Middle Hessian column
    wave_forw=guess; wave_forw(3)=wave_forw(3)+h;
    wave_back=guess; wave_back(3)=wave_back(3)-h;
    [~,~,grad_forw]=grape_xy(wave_forw,spin_system);
    [~,~,grad_back]=grape_xy(wave_back,spin_system);
    grad_forw=squeeze(grad_forw(:,:,1)); grad_forw=grad_forw(:);
    grad_back=squeeze(grad_back(:,:,1)); grad_back=grad_back(:);
    hess_num=(grad_forw-grad_back)/(2*h);
    if norm(hess_anl(:,3)-hess_num,1)<1e-6*norm(hess_num,1)
        disp([formalisms{n} ' middle column test passed']);
    else
        error([formalisms{n} ' middle column test failed']);
    end

end

end
