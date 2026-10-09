% Time-domain Z magnetisation dynamics in the Diels-Alder cycloaddition 
% of acetylene to butadiene, demonstrating the non-linear kinetics module.
% Concentrations occupy unit coordinates; additive product arrival shares
% its unit source equally between reactants. The prescribed concentration
% history and the two-point spin-propagation workflow are retained.
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

% Chemical parts and initial concentrations, mol/L
inter.chem.parts={1:2, 3:8, 9:16, 17:22};
inter.chem.concs=[1e-2 2e-2 0 0.1];

% Additive cycloaddition with rate constant in L/(mol*s)
inter.chem.reactions={struct('reactants',[1 2],'products',3,...
    'matching',[1 9;2 12;3 15;4 16;5 10;6 11;7 14;8 13],...
    'rate',25.0,'closure','additive')};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation={'none', 'none', 'none', 'none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Trace spins for the concentration-only stage of the same reaction network
chem_system=kill_spin(spin_system,1:spin_system.comp.nspins);
K_chem=kinetics(chem_system);

% Time grid (ten seconds)
nsteps=100; tmax=10.0; dt=tmax/nsteps;
time_axis=linspace(0,tmax,nsteps+1); 
 
% Preallocate trajectory 
x=zeros(4,nsteps+1);

% Define initial concentrations
x(:,1)=unit_state(chem_system);
 
% Run Lie group solver
for n=1:nsteps 
    x(:,n+1)=step(chem_system,{@(t,y)1i*K_chem(t,y),n*dt,'LG4'},x(:,n),dt);
end

% Interpolate concentrations as functions of time
A=griddedInterpolant(time_axis,x(1,:),'makima','none');
B=griddedInterpolant(time_axis,x(2,:),'makima','none');
C=griddedInterpolant(time_axis,x(3,:),'makima','none');

% Plot chemical kinetics, excluding ethanol
kfigure(); plot(time_axis',real(x(1:3,:)')); xlim tight; kgrid;
kxlabel('time, seconds'); kylabel('concentration, mol/L');
klegend({'acetylene','butadiene','cyclohexadiene'},...
        'Location','Best');
scale_figure([1.00 0.75]); axis tight; drawnow;

% Compile spin transport and embed prescribed history in unit coordinates
K_spin=kinetics(spin_system);
unit_idx=spin_system.bas.offsets(1:end-1)+1;
unit_embed=sparse(unit_idx,1:4,ones(1,4),spin_system.bas.offsets(end),4);
concs=@(t)[A(t);B(t);C(t);inter.chem.concs(4)];

% Concentration-weighted longitudinal preparation without solvent excitation
eta=state(spin_system,'Lz',1:16);
[~,P]=levelpop('1H',sys.magnet,300);
eta=unit_state(spin_system)+(0.5*P(1)-0.5*P(2))*eta;

% Preallocate the trajectory and get it started
traj=zeros([numel(eta) nsteps+1]); traj(:,1)=eta;

% Run the evolution loop
for n=1:nsteps

    % Keep the user informed
    report(spin_system,['time step ' int2str(n) ...
                        '/' int2str(nsteps)]);

    % Evaluate additive chemistry at the prescribed interval-edge populations
    F_L=1i*K_spin(time_axis(n),unit_embed*concs(time_axis(n)));
    F_R=1i*K_spin(time_axis(n+1),unit_embed*concs(time_axis(n+1)));

    % Take the time step using the two-point Lie quadrature
    traj(:,n+1)=step(spin_system,{F_L,F_R},traj(:,n),dt);

end

% Look at spins in reactants and product
coil=coil_state(spin_system,{'Lz'},{1},'exact');
kfigure(); plot(time_axis,real(coil'*traj)); 
coil=coil_state(spin_system,{'Lz'},{3},'exact');
hold on;  plot(time_axis,real(coil'*traj)); 
coil=coil_state(spin_system,{'Lz'},{9},'exact');
hold on;  plot(time_axis,real(coil'*traj));
xlim tight; kgrid; kxlabel('time, seconds'); 
kylabel('conc.-weighted exp. value, 300K');
klegend({'Acetylene $\hat L_{\rm{Z}}$',...
         'Butadiene $\hat L_{\rm{Z}}$',...
         'Cyclohexadiene $\hat L_{\rm{Z}}$'},...
         'Location','Best');
scale_figure([1.00 0.75]); axis tight;

end

