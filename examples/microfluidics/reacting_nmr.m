% Non-linear reaction kinetics in combination with spin evolution
% (repeated pulse-acquire NMR) and relaxation (Redfield theory).
%
% Calculation time: hours, much faster on GPU.
%
% a.acharya@soton.ac.uk
% madhukar.said@ugent.be
% bruno.linclau@ugent.be
% ilya.kuprov@weizmann.ac.il

function reacting_nmr()

% Import Diels-Alder cycloaddition
[sys,inter,bas,kin]=dac_reaction();

% Magnet field
sys.magnet=14.1;

% Greedy parallelisation
sys.enable={'greedy','gpu'};

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
chem_nsteps=200; chem_tmax=20; 
chem_dt=chem_tmax/chem_nsteps;
chem_time_grid=linspace(0,chem_tmax,chem_nsteps+1); 

% Preallocate concentration trajectory
chem_traj=zeros(5,chem_nsteps+1);

% Initial concentrations, mol/L
chem_traj(:,1)=[0.6; 0.5; 0.0; 0.0; 18.1]; 

% Stage 1: concentration dynamics
for n=1:chem_nsteps 
    chem_traj(:,n+1)=iserstep(spin_system,K,...
                              chem_traj(:,n),(n-1)*chem_dt,chem_dt,'LG4'); 
end

% Plot concentrations, excluding solvent
figure(); plot(chem_time_grid,real(chem_traj(1:4,:))); 
xlim tight; ylim padded; kgrid;
kxlabel('time, seconds'); kylabel('concentration, mol/L');
klegend({'cyclopentadiene','acrylonitrile', ...
         'endo-norbornene carbonitrile',...
         'exo-norbornene carbonitrile'}, ...
         'Location','northeast'); drawnow;

% Interpolate concentrations as functions of time
A=griddedInterpolant(chem_time_grid,chem_traj(1,:),'makima','none');
B=griddedInterpolant(chem_time_grid,chem_traj(2,:),'makima','none');
C=griddedInterpolant(chem_time_grid,chem_traj(3,:),'makima','none');
D=griddedInterpolant(chem_time_grid,chem_traj(4,:),'makima','none');

% Build chemical reaction generators
G1=react_gen(spin_system,kin{1});
G2=react_gen(spin_system,kin{2});

% Get concentration-weighted initial condition, no solvent
eta= A(0)*state(spin_system,'Lz',spin_system.chem.parts{1}) ...
    +B(0)*state(spin_system,'Lz',spin_system.chem.parts{2}) ...
    +C(0)*state(spin_system,'Lz',spin_system.chem.parts{3}) ...
    +D(0)*state(spin_system,'Lz',spin_system.chem.parts{4});
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
    F_L=1i*k1*G1{1}*B(chem_time_grid(n)) ...   % Reaction 1 from substance A
       +1i*k1*A(chem_time_grid(n))*G1{2} ...   % Reaction 1 from substance B
       +1i*k2*G2{1}*B(chem_time_grid(n)) ...   % Reaction 2 from substance A
       +1i*k2*A(chem_time_grid(n))*G2{2};      % Reaction 2 from substance B

    % Build the right interval edge composite evolution generator
    F_R=1i*k1*G1{1}*B(chem_time_grid(n+1)) ... % Reaction 1 from substance A
       +1i*k1*A(chem_time_grid(n+1))*G1{2} ... % Reaction 1 from substance B
       +1i*k2*G2{1}*B(chem_time_grid(n+1)) ... % Reaction 2 from substance A
       +1i*k2*A(chem_time_grid(n+1))*G2{2};    % Reaction 2 from substance B

    % Take the time step using the two-point Lie quadrature
    chem_traj(:,n+1)=step(spin_system,{F_L,F_R},chem_traj(:,n),chem_dt);

end

% Acquisition parameters
parameters.spins={'1H'};
parameters.offset=2328;
parameters.sweep=3500;
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
fids=cell(19,1);

% Acquisitions every second
parfor n=0:18 %#ok<*PFBNS>

    % Pull the initial condition
    eta=chem_traj(:,chem_time_grid==n); 

    % Apply the excitation pulse
    eta=step(spin_system,Hy,eta,pi/2);

    % Get the timing grid
    timing_grid=linspace(n,n+parameters.nsteps*nmr_dt,...
                         parameters.nsteps+1);

    % Get everything to the GPU
    L=gpuArray(H+1i*R);  G11=gpuArray(G1{1});
    G12=gpuArray(G1{2}); G21=gpuArray(G2{1});
    G22=gpuArray(G2{2}); eta=gpuArray(eta); coil=gpuArray(Hp);

    % Get the fid started
    current_fid=gpuArray.zeros(1,parameters.nsteps+1);
    current_fid(1)=hdot(coil,eta);

    % Stage 2: nuclear spin dynamics
    for k=1:parameters.nsteps

        % Keep the user informed
        report(spin_system,['NMR time step ' int2str(k) ...
                            '/' int2str(parameters.nsteps)]);

        % Build the left interval edge composite evolution generator
        F_L=L+1i*k1*G11*B(timing_grid(k)) ...   % Reaction 1 from substance A
             +1i*k1*A(timing_grid(k))*G12 ...   % Reaction 1 from substance B
             +1i*k2*G21*B(timing_grid(k)) ...   % Reaction 2 from substance A
             +1i*k2*A(timing_grid(k))*G22;      % Reaction 2 from substance B

        % Build the right interval edge composite evolution generator
        F_R=L+1i*k1*G11*B(timing_grid(k+1)) ... % Reaction 1 from substance A
             +1i*k1*A(timing_grid(k+1))*G12 ... % Reaction 1 from substance B
             +1i*k2*G21*B(timing_grid(k+1)) ... % Reaction 2 from substance A
             +1i*k2*A(timing_grid(k+1))*G22;    % Reaction 2 from substance B

        % Take the time step using the two-point Lie quadrature
        eta=step(spin_system,{F_L,F_R},eta,nmr_dt);

        % Read out the observable
        current_fid(k+1)=hdot(coil,eta);

    end

    % Store the FID
    fids{n+1}=gather(current_fid);

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

