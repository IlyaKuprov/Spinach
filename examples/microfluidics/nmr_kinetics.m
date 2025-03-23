% Non-linear reaction kinetics in combination with spin evolution
% (repeated pulse-acquire NMR) and relaxation (Redfield theory).
%
% Calculation time: hours.
%
% a.acharya@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function nmr_kinetics()

% Import Diels-Alder cycloaddition
[sys,inter,bas]=dac_reaction();

% Magnet field
sys.magnet=14.1;

% Greedy parallelisation
sys.enable={'greedy'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Rate constants, mol/(L*s)
k1=0.5;  % towards exo  
k2=0.1;  % towards endo

% Cycloaddition reaction generator, including solvent
K=@(t,x)(1i*[-k1*x(2)-k2*x(2)  0                0      0     0;      
              0               -k1*x(1)-k2*x(1)  0      0     0;           
              0                k1*x(1)          0      0     0;
              0                k2*x(1)          0      0     0; 
              0                0                0      0     0]);        

% Kinetic time grid, 20 seconds
nsteps=200; tmax=20; dt=tmax/nsteps;
time_axis=linspace(0,tmax,nsteps+1); 

% Preallocate concentration trajectory
conc_traj=zeros(5,nsteps+1);

% Initial concentrations, mol/L
conc_traj(:,1)=[0.6; 0.5; 0.0; 0.0; 18.1]; 

% Stage 1: concentration dynamics
for n=1:nsteps 
    conc_traj(:,n+1)=iserstep(spin_system,K,conc_traj(:,n),n*dt,dt,'LG4'); 
end

% Plot concentrations, excluding solvent
figure(); plot(time_axis,real(conc_traj(1:4,:))); 
xlim tight; ylim padded; kgrid;
kxlabel('time, seconds'); kylabel('concentration, mol/L');
klegend({'cyclopentadiene','acrylonitrile', ...
         'endo-norbornene carbonitrile',...
         'exo-norbornene carbonitrile'}, ...
         'Location','northeast'); drawnow;

% Interpolate concentrations as functions of time
A=griddedInterpolant(time_axis,conc_traj(1,:),'makima','none');
B=griddedInterpolant(time_axis,conc_traj(2,:),'makima','none');
C=griddedInterpolant(time_axis,conc_traj(3,:),'makima','none');
D=griddedInterpolant(time_axis,conc_traj(4,:),'makima','none');

% Nonlinear kinetics generator
reaction{1}.reactants=[1 2];  % cyclopentadiene and acrylonitrile
reaction{1}.products=3;       % into endo-norbornene carbonitrile
reaction{1}.matching=[1 12; 2 17; 3 18; 4 16; 5 10; 6 11; 7 14; 8 15; 9 13];
reaction{2}.reactants=[1 2];  % cyclopentadiene and acrylonitrile
reaction{2}.products=4;       % into endo-norbornene carbonitrile
reaction{2}.matching=[1 21; 2 26; 3 27; 4 25; 5 19; 6 20; 7 23; 8 24; 9 22]; 

% Build reaction generators and send to GPU
G1=react_gen(spin_system,reaction{1});
G2=react_gen(spin_system,reaction{2});

% Get concentration-weighted initial condition, no solvent
eta= A(0)*state(spin_system,'Lz',spin_system.chem.parts{1}) ...
    +B(0)*state(spin_system,'Lz',spin_system.chem.parts{2}) ...
    +C(0)*state(spin_system,'Lz',spin_system.chem.parts{3}) ...
    +D(0)*state(spin_system,'Lz',spin_system.chem.parts{4});
[~,P]=levelpop('1H',sys.magnet,300);
eta=(0.5*P(1)-0.5*P(2))*eta;

% Run chemistry for 20 seconds
chem_dt=0.1; chem_nsteps=200;

% Preallocate the trajectory and get it started
chem_traj=zeros([numel(eta) chem_nsteps+1]); chem_traj(:,1)=eta;

% Run chemistry
for n=1:chem_nsteps

    % Keep the user informed
    report(spin_system,['chemistry time step ' int2str(n) ...
                        '/' int2str(chem_nsteps)]);

    % Build the left interval edge composite evolution generator
    F_L=1i*k1*G1{1}*B(time_axis(n)) ...   % Reaction 1 from substance A
       +1i*k1*A(time_axis(n))*G1{2} ...   % Reaction 1 from substance B
       +1i*k2*G2{1}*B(time_axis(n)) ...   % Reaction 2 from substance A
       +1i*k2*A(time_axis(n))*G2{2};      % Reaction 2 from substance B

    % Build the right interval edge composite evolution generator
    F_R=1i*k1*G1{1}*B(time_axis(n+1)) ... % Reaction 1 from substance A
       +1i*k1*A(time_axis(n+1))*G1{2} ... % Reaction 1 from substance B
       +1i*k2*G2{1}*B(time_axis(n+1)) ... % Reaction 2 from substance A
       +1i*k2*A(time_axis(n+1))*G2{2};    % Reaction 2 from substance B

    % Take the time step using the two-point Lie quadrature
    chem_traj(:,n+1)=step(spin_system,{F_L,F_R},chem_traj(:,n),chem_dt);

end

% Acquisition parameters
parameters.spins={'1H'};
parameters.offset=2328;
parameters.sweep=3500;
parameters.nsteps=4096;

% Get spin evolution generators 
H=hamiltonian(assume(spin_system,'nmr'));
H=frqoffset(spin_system,H,parameters);
R=relaxation(spin_system);

% Detect transverse magnetisation
coil=state(spin_system,'L+','1H');

% Run NMR every second
chem_traj=chem_traj(:,1:10:end);
start_times=0:(numel(chem_traj)-1);

% Apply the excitation pulse
Hy=operator(spin_system,'Ly','1H');
chem_traj=step(spin_system,Hy,chem_traj,pi/2);

% Timing grid of NMR stage
nmr_dt=1/parameters.sweep;

% Preallocate FID array
fids=cell(numel(start_times),1);

% Run NMR experiments
parfor n=1:size(chem_traj,2) %#ok<*PFBNS>

    % Get the timing grid
    timing_grid=linspace(start_times(n),...
                         start_times(n)+parameters.nsteps*nmr_dt,...
                         parameters.nsteps+1);

    % Gab the initial condition
    eta=chem_traj(:,n);

    % Get the fid started
    current_fid=zeros(1,parameters.nsteps+1);
    current_fid(1)=hdot(coil,eta);

    % Stage 2: nuclear spin dynamics
    for k=1:parameters.nsteps

        % Keep the user informed
        report(spin_system,['NMR time step ' int2str(k) ...
                            '/' int2str(parameters.nsteps)]);

        % Build the left interval edge composite evolution generator
        F_L=H+1i*R+1i*k1*G1{1}*B(timing_grid(k)) ...   % Reaction 1 from substance A
                  +1i*k1*A(timing_grid(k))*G1{2} ...   % Reaction 1 from substance B
                  +1i*k2*G2{1}*B(timing_grid(k)) ...   % Reaction 2 from substance A
                  +1i*k2*A(timing_grid(k))*G2{2};      % Reaction 2 from substance B

        % Build the right interval edge composite evolution generator
        F_R=H+1i*R+1i*k1*G1{1}*B(timing_grid(k+1)) ... % Reaction 1 from substance A
                  +1i*k1*A(timing_grid(k+1))*G1{2} ... % Reaction 1 from substance B
                  +1i*k2*G2{1}*B(timing_grid(k+1)) ... % Reaction 2 from substance A
                  +1i*k2*A(timing_grid(k+1))*G2{2};    % Reaction 2 from substance B

        % Take the time step using the two-point Lie quadrature
        eta=step(spin_system,{F_L,F_R},eta,nmr_dt);

        % Read out the observable
        current_fid(k+1)=hdot(coil,eta);

    end

    % Store the FID
    fids{n}=current_fid;

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

% Export the figure
exportgraphics(gcf,'nmr_kinetics.png','Resolution',1200);
savefig(gcf,'nmr_kinetics.fig');

end

