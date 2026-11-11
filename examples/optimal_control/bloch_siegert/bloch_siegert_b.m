% Bloch-Siegert shift compensation functionality demo. The
% script optimises a universal rotation pulse for a range 
% of resonance offsets. As the control power is increased,
% Bloch-Siegert shift starts to reduce the fidelity unless
% it is correctly accounted for.
%
% Calculation time: many hours (spaghetti plots)
%
% aditya.dev@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il

function bloch_siegert_b()

% Magnet field
sys.magnet=28.18;

% Set isotopes
sys.isotopes={'13C'};

% Set interactions 
inter.zeeman.scalar={0.0};

% Set formalism and basis
bas.formalism='sphten-liouv';
bas.approximation={'none'};

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up spin states
Sx=state(spin_system,'Lx','13C'); Sx=Sx/norm(full(Sx),2);
Sy=state(spin_system,'Ly','13C'); Sy=Sy/norm(full(Sy),2);
Sz=state(spin_system,'Lz','13C'); Sz=Sz/norm(full(Sz),2);

% Get the control operators
Lx=operator(spin_system,'Lx','13C');
Ly=operator(spin_system,'Ly','13C');

% Get the offset operator
Lz=operator(spin_system,'Lz','13C');

% Get the drift Hamiltonian
D=hamiltonian(assume(spin_system,'nmr'));

% Optimal control settings
control.isotopes={'13C'};
control.channels=[1; 1];
control.drifts={{D}};
control.operators={Lx,Ly};
control.off_ops={Lz};
control.offsets={linspace(-3e4,3e4,101)};
control.rho_init={ Sx Sy Sz};
control.rho_targ={-Sz Sy Sx};
control.method='rbfgs';
control.max_iter=500;
control.tol_x=1e-4;
control.plotting={};

% Power levels to sweep (rad/s)
zeeman_frq=abs(spin('13C')*sys.magnet);
pwr_list=zeeman_frq*linspace(1e-3,0.3,50);

% Compute relative control powers
relative_power=pwr_list./zeeman_frq;

% Get a figure going
kfigure(); hold on;
scale_figure([1.0 0.65]);

% Spaghetti plot loop
for k=1:20

    % Preallocate results
    fid_a=zeros(size(pwr_list));
    fid_b=zeros(size(pwr_list));

    % Common initial guess
    guess=randn(2,50)/10;

    % Loop over control power
    for n=1:numel(pwr_list)

        % Set current power
        control.pwr_levels=pwr_list(n);

        % Set pulse slice duration
        dt=(pi/100)/control.pwr_levels;
        control.pulse_dt=dt*ones(1,50);

        % BSS settings: off, on
        control.bsiegert=false();
        setting_a=optimcon(spin_system,control);
        control.bsiegert=true();
        setting_b=optimcon(spin_system,control);

        % Optimisation with and without BSS
        pulse_a=fmaxnewton(setting_a,@grape_xy,guess);
        pulse_b=fmaxnewton(setting_b,@grape_xy,guess);

        % Evaluation for a system with BSS
        [~,fid_a(n)]=ensemble(pulse_a,setting_b);
        [~,fid_b(n)]=ensemble(pulse_b,setting_b);

    end

    % Do the plotting
    plot(relative_power,1-fid_a,'Color',[0.000 0.447 0.741]);
    plot(relative_power,1-fid_b,'Color',[0.850 0.325 0.098]);
    drawnow;

end

% Embellishments
set(gca,'Yscale','log'); xlim tight; ylim padded; 
kxlabel('relative control power $|\omega_1 / \omega_0|$');
klegend({'GRAPE','GRAPE + BSS'},'Location','East');
kylabel('infidelity'); kbox; kletter('B'); kgrid;

% Save the figure
savefig(gcf,'bloch_siegert_b.fig');
pdf_and_png('bloch_siegert_b');

end

