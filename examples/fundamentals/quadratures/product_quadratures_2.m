% A test of Lie-group product quadratures on a chirped frequency 
% oscillator with radiation damping that has a state-dependent
% and time-dependent evolution generator.
%
% a.acharya@soton.ac.uk
% i.kuprov@soton.ac.uk

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

% Run reference RKMK4 simulation and keep the trajectory
np_ref=2^16; dt=0.5/(np_ref-1); mu_traj=zeros(3,np_ref);
mu_traj(:,1)=euler2dcm(0,pi*178/180,0)*[0 0 1]';
for n=2:np_ref 
    mu_traj(:,n)=iserstep(spin_system,G,mu_traj(:,n-1),(n-2)*dt,dt,'RKMK4');
end
mu_ref=mu_traj(:,end);

% Benchmark arrays
np=ceil(2.^linspace(9,15,20)); bench=zeros(numel(np),3);

% Benchmarking loop
for k=1:numel(np)

    % Half a second
    dt=0.5/(np(k)-1);

    % Piecewise-constant, left edge
    mu=euler2dcm(0,pi*178/180,0)*[0 0 1]';
    for n=2:np(k) 
        mu=iserstep(spin_system,G,mu,(n-2)*dt,dt,'PWCL');
    end
    bench(k,1)=norm(mu-mu_ref)/norm(mu_ref);

    % Lie group, order 2
    mu=euler2dcm(0,pi*178/180,0)*[0 0 1]';
    for n=2:np(k) 
        mu=iserstep(spin_system,G,mu,(n-2)*dt,dt,'LG2');
    end
    bench(k,2)=norm(mu-mu_ref)/norm(mu_ref);

    % Lie group, order 4
    mu=euler2dcm(0,pi*178/180,0)*[0 0 1]';
    for n=2:np(k) 
        mu=iserstep(spin_system,G,mu,(n-2)*dt,dt,'LG4');
    end
    bench(k,3)=norm(mu-mu_ref)/norm(mu_ref);

end

% Plotting
figure(); scale_figure([1.5 0.6]);
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
klegend({'LP','LG-2','LG-4'},...
        'Location','SouthWest');
set(gca,'YTick',10.^(-7:2:0));
set(gca,'MinorGridLineStyle','-'); 
set(gca,'MinorGridColor',0.9*[1 1 1]);
set(gca,'GridColor',0.9*[1 1 1]);

end

