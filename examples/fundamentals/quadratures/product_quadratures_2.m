% A test of Lie-group product quadratures on a chirped frequency 
% oscillator with radiation damping that has a state-dependent
% and time-dependent evolution generator.
%
% a.acharya@soton.ac.uk
% ilya.kuprov@weizmann.ac.il 

function product_quadratures_2()

% Set system parameters
cr=2*pi*400;   % 400 Hz/s chirp rate
r1=10; r2=10;  % 10 Hz relaxation
rrd=40;        % 40 Hz radiation damping

% Bootstrap the object
spin_system=bootstrap();

% Make Bloch-Maxwell generator
G=@(t,mu)(-1i*[r2 cr*t 0; -cr*t r2 0; 0 0 r1]...
          -1i*rrd*[mu(3) 0 0; 0 mu(3) 0; mu(1) mu(2) 0]);

% Set initial magnetisation
mu0=euler2dcm(0,pi*178/180,0)*[0 0 1]';

% Run reference RKMK-DP8 simulation and keep the trajectory
np_ref=2^12; dt_ref=0.5/(np_ref-1); mu_traj=zeros(3,np_ref);
mu_traj(:,1)=mu0;
for n=2:np_ref
    t_curr=(n-2)*dt_ref;
    mu_traj(:,n)=iserstep(spin_system,{G,t_curr,'RKMK-DP8'},mu_traj(:,n-1),dt_ref);
end
mu_ref=mu_traj(:,end);

% Benchmark arrays
np=ceil(2.^linspace(8,11,10)); bench=zeros(numel(np),7);

% Benchmarking loop
for k=1:numel(np)

    % Half a second
    dt=0.5/(np(k)-1);

    % Piecewise-constant, left edge
    mu=mu0;
    for n=2:np(k)
        t_curr=(n-2)*dt;
        mu=iserstep(spin_system,{G,t_curr,'PWCL'},mu,dt);
    end
    bench(k,1)=norm(mu-mu_ref)/norm(mu_ref);

    % Lie group, order 2
    mu=mu0;
    for n=2:np(k)
        t_curr=(n-2)*dt;
        mu=iserstep(spin_system,{G,t_curr,'LG2'},mu,dt);
    end
    bench(k,2)=norm(mu-mu_ref)/norm(mu_ref);

    % Lie group, order 4
    mu=mu0;
    for n=2:np(k)
        t_curr=(n-2)*dt;
        mu=iserstep(spin_system,{G,t_curr,'LG4'},mu,dt);
    end
    bench(k,3)=norm(mu-mu_ref)/norm(mu_ref);

    % Lie group, order 4 (from Casas-Iserles Appendix A1)
    mu=mu0;
    for n=2:np(k)
        t_curr=(n-2)*dt;
        mu=iserstep(spin_system,{G,t_curr,'LG4A'},mu,dt);
    end
    bench(k,4)=norm(mu-mu_ref)/norm(mu_ref);

    % Fourth-order RKMK
    mu=mu0;
    for n=2:np(k)
        t_curr=(n-2)*dt;
        mu=iserstep(spin_system,{G,t_curr,'RKMK4'},mu,dt);
    end
    bench(k,5)=norm(mu-mu_ref)/norm(mu_ref);

    % Fifth-order Dormand-Prince RKMK
    mu=mu0;
    for n=2:np(k)
        t_curr=(n-2)*dt;
        mu=iserstep(spin_system,{G,t_curr,'RKMK-DP5'},mu,dt);
    end
    bench(k,6)=norm(mu-mu_ref)/norm(mu_ref);

    % Eighth-order Dormand-Prince RKMK
    mu=mu0;
    for n=2:np(k)
        t_curr=(n-2)*dt;
        mu=iserstep(spin_system,{G,t_curr,'RKMK-DP8'},mu,dt);
    end
    bench(k,7)=norm(mu-mu_ref)/norm(mu_ref);

end

% Plotting
kfigure(); scale_figure([1.8 0.6]);
time_axis=linspace(0,0.5,np_ref);
subplot(1,2,1); plot(time_axis,mu_traj);
kgrid; xlim tight; ylim padded;
kxlabel('time, seconds');
kylabel('magnetisation, a.u.');
klegend({'x-magnetisation','y-magnetisation',...
         'z-magnetisation'},'Location','SouthEast');
subplot(1,2,2); plot(np',bench);
kgrid; xlim tight; ylim([1e-7 1e0]);
set(gca,'YScale','log','XScale','log');
kxlabel('number of points in the time grid');
kylabel('$\|$difference$\|$/$\|$exact$\|$');
klegend({'LP','LG-2','LG-4','LG-4A',...
         'RKMK4','RKMK-DP5','RKMK-DP8'},...
        'Location','SouthWest');
set(gca,'YTick',10.^(-7:2:0));
set(gca,'MinorGridLineStyle','-'); 
set(gca,'MinorGridColor',0.9*[1 1 1]);
set(gca,'GridColor',0.9*[1 1 1]);

% Compute empirical convergence orders
idx=floor(2*numel(np)/3):numel(np);
log_np=log(np(idx));
method_names={'PWCL','LG2','LG4','LG4A','RKMK4','RKMK-DP5','RKMK-DP8'};

% Print empirical convergence orders
fprintf('\nEmpirical convergence orders (slope of log-log error curve):\n');
fprintf('------------------------------------------------------------\n');
for n=1:size(bench,2)
    log_err=log(bench(idx,n));
    p=polyfit(log_np,log_err,1);
    q_empirical=-p(1);
    fprintf('%-10s : %6.3f\n',method_names{n},q_empirical);
end
fprintf('------------------------------------------------------------\n');

end

