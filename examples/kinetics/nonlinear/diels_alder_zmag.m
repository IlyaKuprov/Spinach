% Time-domain Z magnetisation dynamics in the Diels-Alder cycloaddition 
% of acetylene to butadiene, demonstrating the non-linear kinetics module.
%
% Calculation time: minutes.
%
% ilya.kuprov@weizmann.ac.il
% a.acharya@soton.ac.uk

function diels_alder_zmag()

% DFT import options
options.min_j=2.0;         % Minimum J-coupling, Hz
options.style='harmonics'; % Plotting style

% Load and display acetylene      (substance A)
props_a=gparse('acetylene.out');
[sys_a,inter_a]=g2spinach(props_a,{{'H','1H'}},31.8,options);
kfigure(); scale_figure([1.0 1.0]); 
subplot(1,2,1); cst_display(props_a,{'C'},0.005,[],options); 
                camorbit(+45,-45); ktitle('$^{13}$C CST');
subplot(1,2,2); cst_display(props_a,{'H'},0.05,[],options); 
                camorbit(+45,-45); ktitle('$^{1}$H CST'); drawnow;

% Load and display butadiene      (substance B)
props_b=gparse('butadiene.out');
[sys_b,inter_b]=g2spinach(props_b,{{'H','1H'}},31.8,options);
kfigure(); scale_figure([1.5 1.0]);
subplot(1,2,1); cst_display(props_b,{'C'},0.01,[],options); 
                camorbit(+45,-45); ktitle('$^{13}$C CST');
subplot(1,2,2); cst_display(props_b,{'H'},0.1,[],options); 
                camorbit(+45,-45); ktitle('$^{1}$H CST'); drawnow;

% Load and display cyclohexadiene (substance C)
props_c=gparse('cyclohexadiene.out');
[sys_c,inter_c]=g2spinach(props_c,{{'H','1H'}},31.8,options);
kfigure(); scale_figure([1.5 1.0]); 
subplot(1,2,1); cst_display(props_c,{'C'},0.01,[],options); 
                camorbit(+45,-45); ktitle('$^{13}$C CST');
subplot(1,2,2); cst_display(props_c,{'H'},0.1,[],options); 
                camorbit(+45,-45); ktitle('$^{1}$H CST'); drawnow;

% Add natural abundance ethanol   (substance D)
sys_d.isotopes={'1H','1H','1H','1H','1H','1H'};
inter_d.zeeman.matrix={1.26, 1.26, 1.26, ...
                       3.69, 3.69, 2.61}*eye(3);
inter_d.coordinates={[]; []; []; []; []; []};
inter_d.coupling.scalar=zeros(6,6);
inter_d.coupling.scalar(1,[4 5])=7.0;
inter_d.coupling.scalar(2,[4 5])=7.0;
inter_d.coupling.scalar(3,[4 5])=7.0;
inter_d.coupling.scalar=num2cell(inter_d.coupling.scalar);

% Merge the spin systems
[sys,inter]=merge_inp({sys_a,  sys_b,  sys_c,  sys_d},...
                      {inter_a,inter_b,inter_c,inter_d});

% Magnet field
sys.magnet=14.1;

% Chemical parts and unit concentrations
inter.chem.parts={1:2, 3:8, 9:16, 17:22};
inter.chem.concs=[1 1 1 1];

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Initial concentrations, mol/L
A0=1e-2; B0=2e-2; C0=0; D0=0.1;

% 2nd order rate constant
k=25.0;  % mol/(L*s)

% 2nd order reaction generator
K=@(t,x)(1i*[-k*x(2)   0        0        0; 
              0       -k*x(1)   0        0;
              0        k*x(1)   0        0;
              0        0        0        0]);

% Time grid (ten seconds)
nsteps=100; tmax=10.0; dt=tmax/nsteps;
time_axis=linspace(0,tmax,nsteps+1); 
 
% Preallocate trajectory 
x=zeros(4,nsteps+1);

% Define initial concentrations
x(:,1)=[A0 B0 C0 D0]'; 
 
% Run Lie group solver
for n=1:nsteps 
    x(:,n+1)=iserstep(spin_system,K,x(:,n),n*dt,dt,'LG4'); 
end

% Interpolate concentrations as functions of time
A=griddedInterpolant(time_axis,x(1,:),'makima','none');
B=griddedInterpolant(time_axis,x(2,:),'makima','none');
C=griddedInterpolant(time_axis,x(3,:),'makima','none');

% Plot chemical kinetics, excluding ethanol
kfigure(); plot(time_axis',real(x(1:3,:)')); xlim tight; kgrid;
kxlabel('time, seconds'); kylabel('concentration, mol/L');
klegend({'acetylene','butadiene','cyclohexadiene'},...
        'Location','northeast'); 
scale_figure([1.00 0.75]); axis tight; drawnow;

% Build kinetics generators 
reaction.reactants=[1 2];  % which substances are reactants
reaction.products=3;       % which substances are products
reaction.matching=[1  9;  
                   2 12;   % which spin on the left hand
                   3 15;   % side of the reaction arrow
                   4 16;   % goes into which spin on the
                   5 10;   % right hand side
                   6 11;
                   7 14;
                   8 13];
G=react_gen(spin_system,reaction); 

% Get concentration-weighted initial condition, no solvent
eta= A(0)*state(spin_system,'Lz',spin_system.chem.parts{1}) ...
    +B(0)*state(spin_system,'Lz',spin_system.chem.parts{2}) ...
    +C(0)*state(spin_system,'Lz',spin_system.chem.parts{3});
[~,P]=levelpop('1H',sys.magnet,300);
eta=(0.5*P(1)-0.5*P(2))*eta;

% Preallocate the trajectory and get it started
traj=zeros([numel(eta) nsteps+1]); traj(:,1)=eta;

% Run the evolution loop
for n=1:nsteps

    % Keep the user informed
    report(spin_system,['time step ' int2str(n) ...
                        '/' int2str(nsteps)]);

    % Build the left interval edge composite evolution generator
    F_L=1i*k*G{1}*B(time_axis(n))...   % Reaction from substance A
       +1i*k*A(time_axis(n))*G{2};     % Reaction from substance B

    % Build the right interval edge composite evolution generator
    F_R=1i*k*G{1}*B(time_axis(n+1))... % Reaction from substance A
       +1i*k*A(time_axis(n+1))*G{2};   % Reaction from substance B

    % Take the time step using the two-point Lie quadrature
    traj(:,n+1)=step(spin_system,{F_L,F_R},traj(:,n),dt);

end

% Look at spins in reactants and product
coil=state(spin_system,{'Lz'},{1});
kfigure(); plot(time_axis,real(coil'*traj)); 
coil=state(spin_system,{'Lz'},{3});
hold on;  plot(time_axis,real(coil'*traj)); 
coil=state(spin_system,{'Lz'},{9});
hold on;  plot(time_axis,real(coil'*traj));
xlim tight; kgrid; kxlabel('time, seconds'); 
kylabel('conc.-weighted expt. value, 300K');
klegend({'Acetylene $\hat L_{\rm{Z}}$',...
         'Butadiene $\hat L_{\rm{Z}}$',...
         'Cyclohexadiene $\hat L_{\rm{Z}}$'},...
         'Location','northeast'); 
scale_figure([1.00 0.75]); axis tight;

end

