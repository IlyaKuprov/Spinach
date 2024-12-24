% Integrating the Lioville - von Neumann equation through
% one period of the MAS rotor using the piecewise-constant
% Hamiltonian approximation as well as the more accurate
% Lie group integrators.
%
% Calculation time: seconds
%
% a.acharya@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function mas_benchmark()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={5.0,-2.0};
inter.coordinates={[0 0 0]; [0 3.9 0.1]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=+1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system=assume(spin_system,'nmr');

% Get Hamiltonian
[H,Q]=hamiltonian(spin_system);

% Spinning speed, Hz 
rate=50000; T=1/rate;

% Magic angle 
mag=atan(sqrt(2));

% Initial condition
rho_init=state(spin_system,'L+','1H');

% Discretise the period 
np_ref=2^13+1; 
t_ref=linspace(0,T,np_ref);
dt=T/(numel(t_ref)-1);

% Pre-allocate H(t)
HF=cell(1,numel(t_ref));

% Pre-compute Euler angles and H(t)
for n=1:numel(t_ref)
    HF{n}=H+orientation(Q,[0 mag 2*pi*(n-1)/(numel(t_ref)-1)]);

end 

% Run reference simulation
rho_ref=rho_init;
for n=1:2:(numel(t_ref)-2)
    rho_ref=step(spin_system,{HF{n},...
                              HF{n+1},...
                              HF{n+2}},rho_ref,2*dt);
end

% Pre-allocate error array
np=2.^(4:12)+1; error=zeros(numel(np),4);

% Benchmarking loop - maybe parfor 
for k=1:numel(np)

    % Discretise the period
    t=linspace(0,T,np(k));
    dt=T/(numel(t)-1);

    % Pre-allocate H(t)
    HF=cell(1,numel(t));

    % Pre-compute Euler angles and H(t)
    for n=1:numel(t)
        HF{n}=H+orientation(Q,[0 mag 2*pi*(n-1)/(numel(t)-1)]);
    end

    % Piecewise constant, left point 
    rho_one=rho_init;
    for n=1:2:(numel(t)-1)     
        rho_one=step(spin_system,HF{n},rho_one,2*dt);
    end
    error(k,1)=norm(rho_ref-rho_one)/norm(rho_ref);

    % Piecewise constant, mid point 
    rho_one=rho_init;
    for n=1:2:(numel(t)-1)     
        rho_one=step(spin_system,HF{n+1},rho_one,2*dt);
    end
    error(k,2)=norm(rho_ref-rho_one)/norm(rho_ref);

    % Lie group, 2nd order
    rho_two=rho_init;
    for n=1:2:(numel(t)-2)  
        rho_two=step(spin_system,{HF{n}, ...
                                  HF{n+2}},rho_two,2*dt);
    end
    error(k,3)=norm(rho_ref-rho_two)/norm(rho_ref);

    % Lie group, 4th order
    rho_thr=rho_init;
    for n=1:2:(numel(t)-2)  
        rho_thr=step(spin_system,{HF{n}  ...
                                  HF{n+1}...
                                  HF{n+2}},rho_thr,2*dt);
    end
    error(k,4)=norm(rho_ref-rho_thr)/norm(rho_ref);

end

% Plotting
figure(); grid on; plot(np',error);
set(gca,'YScale','log','XScale','log');
set(gca,'MinorGridLineStyle','-'); 
set(gca,'MinorGridColor',0.9*[1 1 1]);
set(gca,'GridColor',0.9*[1 1 1]);
kxlabel('number of points in the rotor grid');
kylabel('$\|$difference$\|$/$\|$exact$\|$'); 
klegend({'LP','MP','LG-2','LG-4'},...
        'Location','SouthWest');
kgrid; xlim tight; ylim([1e-13 1e-2]);
scale_figure([1.00 0.65]);

end


