% Complete microfluidic simulation: diffusion, flow, two second-
% order chemical reactions, and NMR detection in a narrow strip
% of the chip where the coil is assumed to be located. Solvent is a
% spin-bearing pool. The concentration history and spin evolution both
% use the kernel reaction records; the two-stage history workflow is
% retained, with the original frozen-rate concentration time steps.
% Product unit arrival is now shared equally between reactants; this
% preserves the mass-action derivative but changes finite frozen steps.
% It is not numerically identical to the old concentration history.
%
% Calculation time: days, much faster on GPU.
%
% a.acharya@soton.ac.uk
% sylwia.ostrowska@kit.edu
% madhukar.said@ugent.be
% marcel.utz@kit.edu
% bruno.linclau@ugent.be
% ilya.kuprov@weizmann.ac.il

function reacting_flow_nmr()

% Import Diels-Alder cycloaddition
[sys,inter,bas]=dac_reaction();

% Bimolecular rates in L/(mol*s), endo then exo
inter.chem.reactions{1}.rate=2.0;
inter.chem.reactions{2}.rate=1.0;

% Import hydrodynamics information
comsol.mesh_file='chip_mesh.txt';
comsol.velo_file='chip_velo.txt';
comsol.crop={[286.8 287.5],[576.0 579.0]};
comsol.inactivate=[9 10 19 30 20 25 14 13   ...
                   3372 3373 3380 3381 3382 ...
                   3386 3169 3185 3201 3054 ...
                   3077 3055 3053 3078 3186 ...
                   3168 875 899 897 877 876 ...
                   860 858 885 859 883];
mesh=comsol_import(comsol);

% Magnet field
sys.magnet=14.1;

% This needs a GPU
sys.enable={'greedy'}; % 'gpu'

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system.mesh=mesh;

%% Concentration dynamics stage

% Trace all spins for the concentration-only stage, retaining five pools
chem_system=kill_spin(spin_system,1:spin_system.comp.nspins);
K_chem=kinetics(chem_system);

% Strong diffusion
parameters.diff=1e-7;

% Timing parameters
chem_dt=20; chem_nsteps=501;

% Get diffusion and flow generator
GF=flow_gen(spin_system,parameters);

% Concentration trajectory preallocation and the initial state
chem_traj=zeros(5,spin_system.mesh.vor.ncells,chem_nsteps+1);
chem_traj(1,1240,1)=0.50; chem_traj(2,1246,1)=0.25;

% Time evolution loop
for n=1:chem_nsteps

    % Keep the user informed
    report(spin_system,['chemistry + hydrodynamics time step ' int2str(n) ...
                        '/' int2str(chem_nsteps)]);

    % Assemble local chemistry and spatial transport at the current state
    c_curr=chem_traj(:,:,n); c_curr=c_curr(:);
    G=1i*K_chem((n-1)*chem_dt,c_curr)+1i*kron(GF,speye(5));

    % Retain the frozen-rate concentration step of the original workflow
    c_next=step(chem_system,G,c_curr,chem_dt);
    chem_traj(:,:,n+1)=chem_concs(chem_system,c_next).';

end

% Get chemistry time grid
chem_time_grid=linspace(0,chem_dt*chem_nsteps,chem_nsteps+1);

% Concentration functions for each cell
A=cell(spin_system.mesh.vor.ncells,1); B=cell(spin_system.mesh.vor.ncells,1);
C=cell(spin_system.mesh.vor.ncells,1); D=cell(spin_system.mesh.vor.ncells,1);
parfor n=1:spin_system.mesh.vor.ncells
    A{n}=griddedInterpolant(chem_time_grid,squeeze(chem_traj(1,n,:)),'makima','none');
    B{n}=griddedInterpolant(chem_time_grid,squeeze(chem_traj(2,n,:)),'makima','none');
    C{n}=griddedInterpolant(chem_time_grid,squeeze(chem_traj(3,n,:)),'makima','none');
    D{n}=griddedInterpolant(chem_time_grid,squeeze(chem_traj(4,n,:)),'makima','none');
end

%% Full chemistry + hydrodynamics + spin dynamics stage

% Compile the full spin-transport reaction maps once
K_spin=kinetics(spin_system);
unit_idx=spin_system.bas.offsets(1:end-1)+1;
spin_dim=spin_system.bas.offsets(end);
unit_embed=sparse(unit_idx(1:4),1:4,ones(1,4),spin_dim,4);

% Build RF coil phantom
coil_ph=(spin_system.mesh.x(spin_system.mesh.idx.active)>287.0)&...
        (spin_system.mesh.x(spin_system.mesh.idx.active)<287.3)&...
        (spin_system.mesh.y(spin_system.mesh.idx.active)>577.0)&...
        (spin_system.mesh.y(spin_system.mesh.idx.active)<577.5);
coil_ph=double(coil_ph);

% Build control operators
Ly=operator(spin_system,'Ly','1H'); dim=numel(coil_ph);
Ly=polyadic({{spdiags(coil_ph,0,dim,dim),Ly}});

% Build detection states
coil=coil_state(spin_system,'L+','1H','exact');
coil=kron(coil_ph,coil);

% NMR simulation parameters
parameters.spins={'1H'};
parameters.offset=2328;
parameters.sweep=3500;
parameters.nsteps=1024;

% Time step of NMR stage
nmr_dt=1/parameters.sweep;

% Get background evolution generators
H=hamiltonian(assume(spin_system,'nmr'));
H=frqoffset(spin_system,H,parameters);
R=relaxation(spin_system);
F=polyadic({{GF,opium(size(H,1),1)}});
H=polyadic({{opium(size(GF,1),1),H}});
R=polyadic({{opium(size(GF,1),1),R}});

% Move to GPU if requested
if ismember('gpu',sys.enable)
    F=gpuArray(F); H=gpuArray(H); R=gpuArray(R);
end

% Build state operators
LzA=coil_state(spin_system,'Lz',spin_system.chem.parts{1},'exact');
LzB=coil_state(spin_system,'Lz',spin_system.chem.parts{2},'exact');
LzC=coil_state(spin_system,'Lz',spin_system.chem.parts{3},'exact');
LzD=coil_state(spin_system,'Lz',spin_system.chem.parts{4},'exact');

% Preallocate fids array
fids=cell(chem_nsteps,1);

% Parfor prep
n_vals=1:25:chem_nsteps;

% Loop over starting points
parfor j=1:numel(n_vals)

    % dereference
    n=n_vals(j);

    % Prepare only reacting-species magnetisation, leaving solvent spins unexcited
    eta=cell(spin_system.mesh.vor.ncells,1);
    start_time=chem_time_grid(n);
    for k=1:spin_system.mesh.vor.ncells
        eta{k}=A{k}(start_time)*LzA+B{k}(start_time)*LzB+...
               C{k}(start_time)*LzC+D{k}(start_time)*LzD;
        eta{k}(unit_idx)=chem_traj(:,k,n);
    end
    eta=cell2mat(eta);

    % Apply the excitation pulse
    eta=step(spin_system,Ly,eta,pi/2);

    % Get the timing grid
    timing_grid=linspace(start_time,...
                         start_time+parameters.nsteps*nmr_dt,...
                         parameters.nsteps+1);

    % Get the fid started
    current_fid=zeros(1,parameters.nsteps+1);
    current_fid(1)=hdot(coil,eta);

    % Stage 2: nuclear spin dynamics
    for k=1:parameters.nsteps

        % Keep the user informed
        report(spin_system,['NMR time step ' int2str(k) ...
                            '/' int2str(parameters.nsteps)]);

        % Supply the prescribed history through voxel unit coordinates
        concs_left=zeros(4,spin_system.mesh.vor.ncells);
        concs_right=zeros(4,spin_system.mesh.vor.ncells);
        for m=1:spin_system.mesh.vor.ncells %#ok<*PFBNS>
            concs_left(:,m)=[A{m}(timing_grid(k));B{m}(timing_grid(k));...
                             C{m}(timing_grid(k));D{m}(timing_grid(k))];
            concs_right(:,m)=[A{m}(timing_grid(k+1));B{m}(timing_grid(k+1));...
                              C{m}(timing_grid(k+1));D{m}(timing_grid(k+1))];
        end
        eta_left=unit_embed*sparse(concs_left);
        eta_right=unit_embed*sparse(concs_right);
        K_L=K_spin(timing_grid(k),eta_left(:));
        K_R=K_spin(timing_grid(k+1),eta_right(:));

        % Assemble left and right evolution generators
        F_L=H+1i*F+1i*R+1i*K_L; F_R=H+1i*F+1i*R+1i*K_R;

        % Take the time step using the two-point Lie quadrature
        eta=step(spin_system,{F_L,F_R},eta,nmr_dt);

        % Read out the observable
        current_fid(k+1)=hdot(coil,eta);

    end

    % Store the FID
    fids{j}=current_fid;

end

% Merge and apodisation
fids(cellfun(@isempty,fids))=[]; fids=cell2mat(fids);
fids=apodisation(spin_system,fids,{{},{'exp',6'}});

% Zerofilling and Fourier transform
specs=fftshift(fft(fids,16384,2),2);

% Spectrum and time axis ticks
parameters.axis_units='ppm';
parameters.zerofill=16384;
spec_ax=axis_1d(spin_system,parameters);
time_ax=chem_time_grid(1:25:chem_nsteps);

% Waterfall plot
[time_ax,spec_ax]=meshgrid(time_ax,spec_ax); kfigure();
waterfall(time_ax',spec_ax',real(specs),'EdgeColor','k');
kylabel('chemical shift, ppm'); box on;
kxlabel('time, seconds'); kgrid; 
kzlabel('intensity, a.u.'); axis tight;
set(gca,'Projection','perspective');

end

