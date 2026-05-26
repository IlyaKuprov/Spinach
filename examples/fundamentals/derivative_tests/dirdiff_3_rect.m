% Directional derivative test for the Cartesian GRAPE
% module, rectangles integrator.
%
% ilya.kuprov@weizmann.ac.il
% u.rasulov@soton.ac.uk

function dirdiff_3_rect()

% Formalisms to test
formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'};

% Loop over formalisms
for n=1:numel(formalisms)

    % Build the derivative-test system
    [spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalisms{n});

    % Define control parameters
    control.drifts={{H}};                           % Drift
    control.operators={Lx Ly};                      % Controls
    control.rho_init={ Sx Sy Sz};                   % Starting states
    control.rho_targ={-Sz Sy Sx};                   % Target states
    control.pwr_levels=2*pi*linspace(50e3,70e3,10); % Power levels
    control.method='lbfgs';                         % Optimisation method
    control.max_iter=1000;                          % Termination condition
    control.plotting={};                            % Plotting options
    control.integrator='rectangle';                 % Integrator

    % Set the interval grid
    control.pulse_dt=12.8e-6*ones(1,5);

    % Spinach housekeeping
    spin_system=optimcon(spin_system,control);

    % Random guess and finite diff increment
    guess=randn(2,5)/3; h=sqrt(eps('double'));

    % Call GRAPE and request analytical gradient
    [~,~,grad_anl]=grape_xy(guess,spin_system);
    grad_anl=squeeze(grad_anl(:,:,1));

    % Left waveform edge
    wave_forw=guess; wave_forw(1)=wave_forw(1)+h;
    wave_back=guess; wave_back(1)=wave_back(1)-h;
    [~,fid_forw]=grape_xy(wave_forw,spin_system);
    [~,fid_back]=grape_xy(wave_back,spin_system);
    grad_num=(fid_forw(1)-fid_back(1))/(2*h);
    if abs(grad_anl(1)-grad_num)/abs(grad_num)<1e-6
        disp([formalisms{n} ' left edge test passed']);
    else
        error([formalisms{n} ' left edge test failed']);
    end

    % Right waveform edge
    wave_forw=guess; wave_forw(end)=wave_forw(end)+h;
    wave_back=guess; wave_back(end)=wave_back(end)-h;
    [~,fid_forw]=grape_xy(wave_forw,spin_system);
    [~,fid_back]=grape_xy(wave_back,spin_system);
    grad_num=(fid_forw(1)-fid_back(1))/(2*h);
    if abs(grad_anl(end)-grad_num)/abs(grad_num)<1e-6
        disp([formalisms{n} ' right edge test passed']);
    else
        error([formalisms{n} ' right edge test failed']);
    end

    % Waveform midpoint
    wave_forw=guess; wave_forw(3)=wave_forw(3)+h;
    wave_back=guess; wave_back(3)=wave_back(3)-h;
    [~,fid_forw]=grape_xy(wave_forw,spin_system);
    [~,fid_back]=grape_xy(wave_back,spin_system);
    grad_num=(fid_forw(1)-fid_back(1))/(2*h);
    if abs(grad_anl(3)-grad_num)/abs(grad_num)<1e-6
        disp([formalisms{n} ' midpoint test passed']);
    else
        error([formalisms{n} ' midpoint test failed']);
    end

end

end
