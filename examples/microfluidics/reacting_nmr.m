% Non-linear reaction kinetics in combination with spin evolution
% (repeated pulse-acquire NMR) and relaxation (Redfield theory).
% The same additive reaction records generate concentration and spin
% transport; product unit arrival is shared equally between reactants.
% Solvent retains its three protons and remains unexcited.
%
% Calculation time: hours, much faster on GPU.
%
% a.acharya@soton.ac.uk
% madhukar.said@ugent.be
% bruno.linclau@ugent.be
% ilya.kuprov@weizmann.ac.il

function reacting_nmr()

% Import Diels-Alder cycloaddition
[sys,inter,bas]=dac_reaction();

% Magnet field
sys.magnet=14.1;

% Greedy parallelisation
sys.enable={'greedy'}; % 'gpu'

% Bimolecular rates in L/(mol*s), endo then exo
inter.chem.reactions{1}.rate=0.5;
inter.chem.reactions{2}.rate=0.1;
inter.chem.concs=[0.6 0.5 0.0 0.0 18.1];

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Trace spins for concentration dynamics with the same reaction network
chem_system=kill_spin(spin_system,1:spin_system.comp.nspins);
K_chem=kinetics(chem_system);

% Kinetic time grid, 20 seconds
chem_nsteps=200; chem_tmax=20; 
chem_dt=chem_tmax/chem_nsteps;
chem_time_grid=linspace(0,chem_tmax,chem_nsteps+1); 

% Preallocate concentration trajectory
chem_traj=zeros(5,chem_nsteps+1);

% Initial concentrations, mol/L
chem_traj(:,1)=unit_state(chem_system);

% Stage 1: concentration dynamics
for n=1:chem_nsteps 
    chem_traj(:,n+1)=step(chem_system,{@(t,y)1i*K_chem(t,y),(n-1)*chem_dt,'LG4'},...
                          chem_traj(:,n),chem_dt); 
end

% Plot concentrations, excluding solvent
kfigure(); plot(chem_time_grid,real(chem_traj(1:4,:))); 
xlim tight; ylim padded; kgrid;
kxlabel('time, seconds'); kylabel('concentration, mol/L');
klegend({'cyclopentadiene','acrylonitrile', ...
         'endo-norbornene carbonitrile',...
         'exo-norbornene carbonitrile'}, ...
         'Location','Best'); drawnow;

% Interpolate concentrations as functions of time
A=griddedInterpolant(chem_time_grid,chem_traj(1,:),'makima','none');
B=griddedInterpolant(chem_time_grid,chem_traj(2,:),'makima','none');
C=griddedInterpolant(chem_time_grid,chem_traj(3,:),'makima','none');
D=griddedInterpolant(chem_time_grid,chem_traj(4,:),'makima','none');

% Compile spin transport and embed prescribed concentrations in unit coordinates
K_spin=kinetics(spin_system);
unit_idx=spin_system.bas.offsets(1:end-1)+1;
unit_embed=sparse(unit_idx,1:5,ones(1,5),spin_system.bas.offsets(end),5);
concs=@(t)[A(t);B(t);C(t);D(t);inter.chem.concs(5)];

% Concentration-weighted longitudinal preparation without solvent excitation
eta=state(spin_system,'Lz',[spin_system.chem.parts{1:4}]);
[~,P]=levelpop('1H',sys.magnet,300);
eta=unit_state(spin_system)+(0.5*P(1)-0.5*P(2))*eta;

% Preallocate the trajectory and get it started
chem_traj=zeros([numel(eta) chem_nsteps+1]); chem_traj(:,1)=eta;

% Run chemistry
for n=1:chem_nsteps

    % Keep the user informed
    report(spin_system,['chemistry time step ' int2str(n) ...
                        '/' int2str(chem_nsteps)]);

    % Evaluate additive chemistry at prescribed interval-edge populations
    F_L=1i*K_spin(chem_time_grid(n),unit_embed*concs(chem_time_grid(n)));
    F_R=1i*K_spin(chem_time_grid(n+1),unit_embed*concs(chem_time_grid(n+1)));

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
Hp=coil_state(spin_system,'L+','1H','exact');

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

    % Move to GPU if requested
    if ismember('gpu',spin_system.sys.enable)
        L=gpuArray(H+1i*R); eta=gpuArray(eta); coil=gpuArray(Hp);
        current_fid=gpuArray.zeros(1,parameters.nsteps+1);
    else
        L=H+1i*R; coil=Hp;
        current_fid=zeros(1,parameters.nsteps+1);
    end

    % Get the fid started
    current_fid(1)=hdot(coil,eta);

    % Stage 2: nuclear spin dynamics
    for k=1:parameters.nsteps

        % Keep the user informed
        report(spin_system,['NMR time step ' int2str(k) ...
                            '/' int2str(parameters.nsteps)]);

        % Assemble interval-edge chemistry from the prescribed unit populations
        K_L=K_spin(timing_grid(k),unit_embed*concs(timing_grid(k)));
        K_R=K_spin(timing_grid(k+1),unit_embed*concs(timing_grid(k+1)));
        if ismember('gpu',spin_system.sys.enable)
            K_L=gpuArray(K_L); K_R=gpuArray(K_R);
        end
        F_L=L+1i*K_L; F_R=L+1i*K_R;

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
[time_ax,spec_ax]=meshgrid(time_ax,spec_ax); kfigure();
waterfall(time_ax',spec_ax',real(specs),'EdgeColor','k');
kylabel('chemical shift, ppm'); box on;
kxlabel('time, seconds'); kgrid; 
kzlabel('intensity, a.u.'); axis tight;
set(gca,'Projection','perspective');

end

