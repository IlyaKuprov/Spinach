% Benchmarks iserstep higher-order methods on a chirped-frequency oscillator
% with radiation damping, that has a state-dependent, and a time-dependent
% evolution generator. Syntax:
%
%               iserstep_bench_hiord()
%
% Outputs:
%
%     (none) - produces a set of diagnostic plots, and prints the empirical
%              convergence orders for each method
%
% a.acharya@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=iserstep_bench_hiord.m>

function iserstep_bench_hiord()

% Set the chirp rate
cr=2*pi*400;

% Set the relaxation rates
r1=10; r2=10;

% Set the radiation damping rate
rrd=40;

% Bootstrap the object
spin_system=bootstrap();

% Make Bloch-Maxwell generator (Liouvillian, including -1i factors)
G=@(t,mu)(-1i*[r2 cr*t 0; -cr*t r2 0; 0 0 r1]...
          -1i*rrd*[mu(3) 0 0; 0 mu(3) 0; mu(1) mu(2) 0]);

% Set the initial magnetisation
mu0=euler2dcm(0,pi*178/180,0)*[0 0 1]';

% Set the reference grid size
np_ref=2^12;

% Set the reference time step
dt_ref=0.5/(np_ref-1);

% Preallocate the reference trajectory array
mu_traj=zeros(3,np_ref);

% Set the initial state
mu_traj(:,1)=mu0;

% Compute the reference trajectory using RKMK-DP8
for n=2:np_ref
    t_curr=(n-2)*dt_ref;
    mu_traj(:,n)=step(spin_system,{G,t_curr,'RKMK-DP8'},mu_traj(:,n-1),dt_ref);
end

% Extract the reference final state
mu_ref=mu_traj(:,end);

% Set the benchmark grid sizes
np=ceil(2.^linspace(8,10.5,10));

% Preallocate the accuracy array
bench=zeros(numel(np),7);

% Benchmarking loop
for k=1:numel(np)

    % Set the time step
    dt=0.5/(np(k)-1);

    % Run PWCL propagation and time it
    mu=run_method(spin_system,G,mu0,dt,np(k),'PWCL',false);
    bench(k,1)=norm(mu-mu_ref)/norm(mu_ref);


    % Run LG2 propagation and time it
    mu=run_method(spin_system,G,mu0,dt,np(k),'LG2',false);
    bench(k,2)=norm(mu-mu_ref)/norm(mu_ref);


    % Run LG4 propagation and time it
    mu=run_method(spin_system,G,mu0,dt,np(k),'LG4',false);
    bench(k,3)=norm(mu-mu_ref)/norm(mu_ref);


    % Run LG4A propagation and time it
    mu=run_method(spin_system,G,mu0,dt,np(k),'LG4A',false);
    bench(k,4)=norm(mu-mu_ref)/norm(mu_ref);


    % Run RKMK4 propagation and time it
    mu=run_method(spin_system,G,mu0,dt,np(k),'RKMK4',true);
    bench(k,5)=norm(mu-mu_ref)/norm(mu_ref);


    % Run RKMK-DP5 propagation and time it
    mu=run_method(spin_system,G,mu0,dt,np(k),'RKMK-DP5',true);
    bench(k,6)=norm(mu-mu_ref)/norm(mu_ref);


    % Run RKMK-DP8 propagation and time it
    mu=run_method(spin_system,G,mu0,dt,np(k),'RKMK-DP8',true);
    bench(k,7)=norm(mu-mu_ref)/norm(mu_ref);

end

% Plot the reference trajectory
kfigure(); scale_figure([1.8 0.6]);

% Make the reference time axis
time_axis=linspace(0,0.5,np_ref);

% Plot the reference trajectory
subplot(1,2,1);
plot(time_axis,mu_traj);
kgrid; xlim tight; ylim padded;
kxlabel('time, seconds');
kylabel('magnetisation, a.u.');
klegend({'x-magnetisation','y-magnetisation',...
         'z-magnetisation'},'Location','SouthEast');

% Plot error against the grid size
subplot(1,2,2);
plot(np',bench);
kgrid; xlim tight; ylim([1e-7 1e0]);
set(gca,'YScale','log','XScale','log');
kxlabel('number of points in the time grid');
kylabel('relative error');
klegend({'LP','LG-2','LG-4','LG-4A',...
         'RKMK4','RKMK-DP5','RKMK-DP8'},...
        'Location','SouthWest');
set(gca,'YTick',10.^(-7:2:0));
set(gca,'MinorGridLineStyle','-');
set(gca,'MinorGridColor',0.9*[1 1 1]);
set(gca,'GridColor',0.9*[1 1 1]);

% Use last third of data points (asymptotic regime)
idx=floor(2*numel(np)/3):numel(np);

% Take logarithms of the grid sizes
log_np=log(np(idx));

% Set the method names for printing
method_names={...
    'PWCL',...
    'LG2',...
    'LG4',...
    'LG4A',...
    'RKMK4',...
    'RKMK-DP5',...
    'RKMK-DP8'};

% Print the convergence order header
fprintf('\nEmpirical convergence orders (slope of log-log error curve):\n');
fprintf('------------------------------------------------------------\n');

% Loop over the methods
for n=1:size(bench,2)

    % Take logarithms of the errors
    log_err=log(bench(idx,n));

    % Fit log(err)=a*log(np)+b
    p=polyfit(log_np,log_err,1);

    % Compute the empirical order
    q_emp=-p(1);

    % Print the results
    fprintf('%-10s : %6.3f\n',method_names{n},q_emp);

end

% Print the convergence order footer
fprintf('------------------------------------------------------------\n');

end

% Runs a propagation using a specific method
function mu=run_method(spin_system,G,mu0,dt,np,method,use_step)

mu=mu0;
for n=2:np
    t_curr=(n-2)*dt;
    if use_step
        mu=step(spin_system,{G,t_curr,method},mu,dt);
    else
        mu=iserstep(spin_system,{G,t_curr,method},mu,dt);
    end
end

end
