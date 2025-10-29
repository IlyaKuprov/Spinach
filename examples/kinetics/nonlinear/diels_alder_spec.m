% Repeated pulse-acquire experiment during the Diels-Alder cyclo-
% addition of acetylene to butadiene, demonstrating the non-linear
% kinetics module.
%
% Calculation time: hours, much faster on GPU.
%
% ilya.kuprov@weizmann.ac.il
% a.acharya@soton.ac.uk

function diels_alder_spec()

% DFT import options
options.min_j=2.0;         % Minimum J-coupling, Hz
options.style='harmonics'; % Plotting style

% Load and display acetylene      (substance A)
props_a=gparse('acetylene.out');
[sys_a,inter_a]=g2spinach(props_a,{{'H','1H'}},31.8,options);
figure(); scale_figure([1.0 1.0]); 
subplot(1,2,1); cst_display(props_a,{'C'},0.005,[],options); 
                camorbit(+45,-45); ktitle('$^{13}$C CST');
subplot(1,2,2); cst_display(props_a,{'H'},0.05,[],options); 
                camorbit(+45,-45); ktitle('$^{1}$H CST'); drawnow;

% Load and display butadiene      (substance B)
props_b=gparse('butadiene.out');
[sys_b,inter_b]=g2spinach(props_b,{{'H','1H'}},31.8,options);
figure(); scale_figure([1.5 1.0]);
subplot(1,2,1); cst_display(props_b,{'C'},0.01,[],options); 
                camorbit(+45,-45); ktitle('$^{13}$C CST');
subplot(1,2,2); cst_display(props_b,{'H'},0.1,[],options); 
                camorbit(+45,-45); ktitle('$^{1}$H CST'); drawnow;

% Load and display cyclohexadiene (substance C)
props_c=gparse('cyclohexadiene.out');
[sys_c,inter_c]=g2spinach(props_c,{{'H','1H'}},31.8,options);
figure(); scale_figure([1.5 1.0]); 
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

% This needs a GPU
sys.enable={'greedy'}; % 'gpu'

% Relaxation theory parameters
inter.relaxation={'redfield','t1_t2'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={ 1e-12 ... % Acetylene
             20e-12 ... % Butadiene
             50e-12 ... % Cyclohexadiene
              5e-12};   % Ethanol
inter.r1_rates=cell(22,1); inter.r1_rates(:)={0};
inter.r2_rates=cell(22,1); inter.r2_rates(:)={0};
inter.r1_rates(17:22)={0.5}; % Solvent
inter.r2_rates(17:22)={0.5}; % Solvent

% Chemical parts and unit concentrations
inter.chem.parts={1:2, 3:8, 9:16, 17:22};
inter.chem.concs=[1 1 1 1];

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% 2nd order rate constant
k=25.0;  % mol/(L*s)

% 2nd order reaction generator
K=@(t,x)(1i*[-k*x(2)   0        0        0; 
              0       -k*x(1)   0        0;
              0        k*x(1)   0        0;
              0        0        0        0]);

% Kinetic time grid, 10 seconds
chem_nsteps=100; chem_tmax=10; 
chem_dt=chem_tmax/chem_nsteps;
chem_time_grid=linspace(0,chem_tmax,chem_nsteps+1); 

% Preallocate concentration trajectory
conc_traj=zeros(4,chem_nsteps+1);

% Initial concentrations, mol/L
conc_traj(:,1)=[1e-2; 2e-2; 0.0; 17.1]; 

% Stage 1: concentration dynamics
for n=1:chem_nsteps 
    conc_traj(:,n+1)=iserstep(spin_system,K,...
                              conc_traj(:,n),(n-1)*chem_dt,chem_dt,'LG4'); 
end

% Plot concentrations, excluding solvent
figure(); plot(chem_time_grid,real(conc_traj(1:3,:))); 
xlim tight; ylim padded; kgrid;
kxlabel('time, seconds'); kylabel('concentration, mol/L');
klegend({'acetylene','butadiene','cyclohexadiene'}, ...
         'Location','northeast'); 
scale_figure([1.00 0.75]); axis tight; drawnow;

% Interpolate concentrations as functions of time
A=griddedInterpolant(chem_time_grid,conc_traj(1,:),'makima','none');
B=griddedInterpolant(chem_time_grid,conc_traj(2,:),'makima','none');
C=griddedInterpolant(chem_time_grid,conc_traj(3,:),'makima','none');

% Build chemical reaction generators
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
chem_traj=zeros([numel(eta) chem_nsteps+1]); chem_traj(:,1)=eta;

% Run chemistry
for n=1:chem_nsteps

    % Keep the user informed
    report(spin_system,['chemistry time step ' int2str(n) ...
                        '/' int2str(chem_nsteps)]);

    % Build the left interval edge composite evolution generator
    F_L=1i*k*B(chem_time_grid(n))*G{1}   ... % Reaction from substance A
       +1i*k*A(chem_time_grid(n))*G{2};      % Reaction from substance B

    % Build the right interval edge composite evolution generator
    F_R=1i*k*B(chem_time_grid(n+1))*G{1} ... % Reaction from substance A
       +1i*k*A(chem_time_grid(n+1))*G{2};    % Reaction from substance B

    % Take the time step using the two-point Lie quadrature
    chem_traj(:,n+1)=step(spin_system,{F_L,F_R},chem_traj(:,n),chem_dt);

end

% Acquisition parameters
parameters.spins={'1H'};
parameters.offset=2370;
parameters.sweep=4000;
parameters.nsteps=4096;

% Time step of NMR stage
nmr_dt=1/parameters.sweep;

% Get spin evolution generators 
H=hamiltonian(assume(spin_system,'nmr'));
H=frqoffset(spin_system,H,parameters);
R=relaxation(spin_system);

% Get the pulse operator
Hy=operator(spin_system,'Ly','1H');

% Detect transverse magnetisation
Hp=state(spin_system,'L+','1H');

% Preallocate FID array
fids=cell(9,1);

% Run NMR experiments
parfor n=0:8 %#ok<*PFBNS>

    % Pull the initial condition
    eta=chem_traj(:,chem_time_grid==n); 

    % Apply the excitation pulse
    eta=step(spin_system,Hy,eta,pi/2);

    % Get the timing grid
    timing_grid=linspace(n,n+parameters.nsteps*nmr_dt,...
                         parameters.nsteps+1);

    % Get everything to the GPU
    L=gpuArray(H+1i*R); G1=gpuArray(G{1});
    G2=gpuArray(G{2}); eta=gpuArray(eta); coil=gpuArray(Hp);

    % Get the fid started
    current_fid=gpuArray.zeros(1,parameters.nsteps+1);
    current_fid(1)=hdot(coil,eta);

    % Stage 2: nuclear spin dynamics
    for k=1:parameters.nsteps

        % Keep the user informed
        report(spin_system,['NMR time step ' int2str(k) ...
                            '/' int2str(parameters.nsteps)]);

        % Build the left interval edge composite evolution generator
        F_L=L+1i*k*G1*B(timing_grid(k))   ... % Reaction from substance A
             +1i*k*A(timing_grid(k))*G2;      % Reaction from substance B

        % Build the right interval edge composite evolution generator
        F_R=L+1i*k*G1*B(timing_grid(k+1)) ... % Reaction from substance A
             +1i*k*A(timing_grid(k+1))*G2;    % Reaction from substance B

        % Take the time step using the two-point Lie quadrature
        eta=step(spin_system,{F_L,F_R},eta,nmr_dt);

        % Read out the observable
        current_fid(k+1)=hdot(coil,eta);

    end

    % Store the FID
    fids{n+1}=current_fid;

end

% Merge and apodisation
fids=cell2mat(fids);
fids=apodisation(spin_system,fids,{{},{'exp',6'}});

% Zerofilling and Fourier transform
specs=fftshift(fft(fids,16384,2),2);

% Spectrum and time axis ticks
parameters.axis_units='ppm';
parameters.zerofill=16384;
spec_ax=axis_1d(spin_system,parameters);
time_ax=(1:size(fids,1))-1;

% Waterfall plot
[time_ax,spec_ax]=meshgrid(time_ax,spec_ax); figure();
waterfall(time_ax',spec_ax',real(specs),'EdgeColor','k');
kylabel('chemical shift, ppm'); box on;
kxlabel('time, seconds'); kgrid; 
kzlabel('intensity, a.u.'); axis tight;
set(gca,'Projection','perspective');

end

