% GRAPE phase Hessian test against finite-differenced gradients,
% rectangles integrator.
%
% ilya.kuprov@weizmann.ac.il
% u.rasulov@soton.ac.uk

function dirdiff_8_rect()

% Formalisms to test
formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'};

% Loop over formalisms
for n=1:numel(formalisms)

    % Build the derivative-test system
    [spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalisms{n});

    % Define control parameters
    control.drifts={{H}};                                      % Drift
    control.operators={Lx,Ly,0.4*Lx+0.2*Ly,0.7*Ly-0.1*Lx};      % Controls
    control.rho_init={ Sx Sy Sz};                              % Starting states
    control.rho_targ={-Sz Sy Sx};                              % Target states
    control.pwr_levels=2*pi*linspace(50e3,70e3,10);            % Power levels
    control.method='newton';                                   % Optimisation method
    control.max_iter=1000;                                     % Termination condition
    control.plotting={};                                       % Plotting options
    control.integrator='rectangle';                            % Integrator

    % Set the interval grid
    control.pulse_dt=12.8e-6*ones(1,5);
    control.amplitudes=[ones(1,5); 0.8+0.1*(1:5)];

    % Spinach housekeeping
    spin_system=optimcon(spin_system,control);

    % Random phases and finite diff increment
    guess=randn(2,5)/3; h=1e-5;

    % Call GRAPE and request analytical Hessian
    [~,~,~,hess_anl]=grape_phase(guess,spin_system);
    hess_anl=squeeze(hess_anl(:,:,1));

    % Leftmost Hessian column
    wave_forw=guess; wave_forw(1)=wave_forw(1)+h;
    wave_back=guess; wave_back(1)=wave_back(1)-h;
    [~,~,grad_forw]=grape_phase(wave_forw,spin_system);
    [~,~,grad_back]=grape_phase(wave_back,spin_system);
    grad_forw=squeeze(grad_forw(:,:,1)); grad_forw=grad_forw(:);
    grad_back=squeeze(grad_back(:,:,1)); grad_back=grad_back(:);
    hess_num=(grad_forw-grad_back)/(2*h);
    if norm(hess_anl(:,1)-hess_num,1)<1e-5*norm(hess_num,1)
        disp([formalisms{n} ' leftmost column test passed']);
    else
        error([formalisms{n} ' leftmost column test failed']);
    end

    % Rightmost Hessian column
    wave_forw=guess; wave_forw(end)=wave_forw(end)+h;
    wave_back=guess; wave_back(end)=wave_back(end)-h;
    [~,~,grad_forw]=grape_phase(wave_forw,spin_system);
    [~,~,grad_back]=grape_phase(wave_back,spin_system);
    grad_forw=squeeze(grad_forw(:,:,1)); grad_forw=grad_forw(:);
    grad_back=squeeze(grad_back(:,:,1)); grad_back=grad_back(:);
    hess_num=(grad_forw-grad_back)/(2*h);
    if norm(hess_anl(:,end)-hess_num,1)<1e-5*norm(hess_num,1)
        disp([formalisms{n} ' rightmost column test passed']);
    else
        error([formalisms{n} ' rightmost column test failed']);
    end

    % Middle Hessian column
    wave_forw=guess; wave_forw(5)=wave_forw(5)+h;
    wave_back=guess; wave_back(5)=wave_back(5)-h;
    [~,~,grad_forw]=grape_phase(wave_forw,spin_system);
    [~,~,grad_back]=grape_phase(wave_back,spin_system);
    grad_forw=squeeze(grad_forw(:,:,1)); grad_forw=grad_forw(:);
    grad_back=squeeze(grad_back(:,:,1)); grad_back=grad_back(:);
    hess_num=(grad_forw-grad_back)/(2*h);
    if norm(hess_anl(:,5)-hess_num,1)<1e-5*norm(hess_num,1)
        disp([formalisms{n} ' middle column test passed']);
    else
        error([formalisms{n} ' middle column test failed']);
    end

end
end
