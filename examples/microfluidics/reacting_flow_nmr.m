% Complete microfuidic simulation: diffusion, flow, two second-
% order chemical reactions, and NMR detection in a narrow strip
% of the chip where the coil is assumed to be located.
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
[sys,inter,bas,kin]=dac_reaction();

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

% Rate constants, mol/(L*s)
k1=2.0;  % towards exo  
k2=1.0;  % towards endo

% Cycloaddition reaction generator, including solvent
K=@(x)([-k1*x(2)-k2*x(2)  0                0      0     0;      
         0               -k1*x(1)-k2*x(1)  0      0     0;           
         0                k1*x(1)          0      0     0;
         0                k2*x(1)          0      0     0; 
         0                0                0      0     0]);  

% Strong diffusion (TODO: needs to be different for each substance)
parameters.diff=1e-7;

% Timing parameters
chem_dt=20; chem_nsteps=501;

% Get diffusion and flow generator
GF=flow_gen(spin_system,parameters);

% Concentratoin trajectory preallocation and the initial state
chem_traj=zeros(5,spin_system.mesh.vor.ncells,chem_nsteps+1);
chem_traj(1,1240,1)=0.50; chem_traj(2,1246,1)=0.25;

% Time evolution loop
for n=1:chem_nsteps

    % Keep the user informed
    report(spin_system,['chemistry + hydrodynamics time step ' int2str(n) ...
                        '/' int2str(chem_nsteps)]);

    % Build a kinetics generator in each cell
    GK=zeros([5 5 spin_system.mesh.vor.ncells],'like',1i);
    parfor k=1:spin_system.mesh.vor.ncells
        GK(:,:,k)=K(chem_traj(:,k,n));
    end

    % Assemble the evolution generator
    G=1i*spblkdiag(GK)+1i*kron(GF,speye(5));

    % Take the time step
    c_curr=chem_traj(:,:,n); c_curr=c_curr(:);
    c_next=step(spin_system,G,c_curr,chem_dt);
    chem_traj(:,:,n+1)=reshape(c_next,[5 spin_system.mesh.vor.ncells]);

end

% Get chemistry time grid
chem_time_grid=linspace(0,chem_dt*chem_nsteps,chem_nsteps+1);

% Concentration functions for each cell
A=cell(spin_system.mesh.vor.ncells,1); B=cell(spin_system.mesh.vor.ncells,1);
C=cell(spin_system.mesh.vor.ncells,1); D=cell(spin_system.mesh.vor.ncells,1);
parfor n=1:spin_system.mesh.vor.ncells
    A{n}=griddedInterpolant(chem_time_grid,squeeze(chem_traj(1,n,:)),'makima','none');
    B{n}=griddedInterpolant(chem_time_grid,squeeze(chem_traj(1,n,:)),'makima','none');
    C{n}=griddedInterpolant(chem_time_grid,squeeze(chem_traj(1,n,:)),'makima','none');
    D{n}=griddedInterpolant(chem_time_grid,squeeze(chem_traj(1,n,:)),'makima','none');
end

%% Full chemistry + hydrodynamics + spin dynamics stage

% Build chemical reaction generators
G1=react_gen(spin_system,kin{1});
G2=react_gen(spin_system,kin{2});

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
coil=state(spin_system,'L+','1H');
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
F=gpuArray(polyadic({{GF,opium(size(H,1),1)}}));
H=gpuArray(polyadic({{opium(size(GF,1),1),H}}));
R=gpuArray(polyadic({{opium(size(GF,1),1),R}}));

% Build state operators
LzA=state(spin_system,'Lz',spin_system.chem.parts{1});
LzB=state(spin_system,'Lz',spin_system.chem.parts{2});
LzC=state(spin_system,'Lz',spin_system.chem.parts{3});
LzD=state(spin_system,'Lz',spin_system.chem.parts{4});

% Preallocate fieds array
fids=cell(chem_nsteps,1);

% Parfor prep
n_vals=1:25:chem_nsteps;

% Loop over starting points
parfor j=1:numel(n_vals)

    % dereference
    n=n_vals(j);

    % Build the initial condition
    eta=cell(spin_system.mesh.vor.ncells,1);
    start_time=chem_time_grid(n);
    for k=1:spin_system.mesh.vor.ncells
        eta{k}=A{k}(start_time)*LzA+B{k}(start_time)*LzB+...
               C{k}(start_time)*LzC+D{k}(start_time)*LzD;    
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

        % Build the kinetics generator
        K_L=cell(spin_system.mesh.vor.ncells,1);
        K_R=cell(spin_system.mesh.vor.ncells,1);
        for m=1:spin_system.mesh.vor.ncells %#ok<*PFBNS>

            % Build the left interval edge kinetics generator
            K_L{m}=k1*G1{1}*B{m}(timing_grid(k))+ ...   % Reaction 1 from substance A
                   k1*A{m}(timing_grid(k))*G1{2}+ ...   % Reaction 1 from substance B
                   k2*G2{1}*B{m}(timing_grid(k))+ ...   % Reaction 2 from substance A
                   k2*A{m}(timing_grid(k))*G2{2};       % Reaction 2 from substance B

            % Build the right interval edge composite evolution generator
            K_R{m}=k1*G1{1}*B{m}(timing_grid(k+1))+ ... % Reaction 1 from substance A
                   k1*A{m}(timing_grid(k+1))*G1{2}+ ... % Reaction 1 from substance B
                   k2*G2{1}*B{m}(timing_grid(k+1))+...  % Reaction 2 from substance A
                   k2*A{m}(timing_grid(k+1))*G2{2};     % Reaction 2 from substance B

        end
        K_L=matlab.internal.math.blkdiag(K_L{:});
        K_R=matlab.internal.math.blkdiag(K_R{:});

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
[time_ax,spec_ax]=meshgrid(time_ax,spec_ax); figure();
waterfall(time_ax',spec_ax',real(specs),'EdgeColor','k');
kylabel('chemical shift, ppm'); box on;
kxlabel('time, seconds'); kgrid; 
kzlabel('intensity, a.u.'); axis tight;
set(gca,'Projection','perspective');

end

