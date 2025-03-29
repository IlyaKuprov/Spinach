% First draft of the magnetohydrodynamics context for microfluidic simu-
% lations. Generates evolution generators and passes them on to the pul-
% se sequence function, which should be supplied as a handle. Syntax:
%
%         answer=meshflow(spin_system,pulse_sequence,parameters)
%
% Parameters:
%
%  pulse_sequence     - pulse sequence function handle. See the
%                       experiments directory for the list of
%                       pulse sequences that ship with Spinach.
%
% The following phantoms must be specified: hamiltonian, relaxation, ki-
% netics, initial condition, detection state. Operator phantoms must be
% specified in the following way:
%
%                   parameters.R_ph={Ph1,Ph2,...,PhN}
%                   parameters.R_op={R1,R2,...,RN}
%
% where PhN have the same dimension as the sample voxel grid and RN are
% relaxation superoperators. Likewise for the following:
%
%                   parameters.K_ph, parameters.K_op
%                   parameters.H_ph, parameters.H_op
%
% The initial condition phantom reflects the fact that different voxels
% might start off in a different spin state. It must be specified in the
% following way:
%
%                 parameters.rho0_ph={Ph1,Ph2,...,PhN}
%                 parameters.rho0_op={rho1,rho2,...,rhoN}
%
% where PhN have the same dimension as the sample voxel grid and rhoN are 
% spin states obtained from state() function.
%
% The detection state phantom reflects the fact that different voxels mi-
% ght be detected at different angles and with different sensitivity. It
% must be specified in the following way:
%
%                 parameters.coil_ph={Ph1,Ph2,...,PhN}
%                 parameters.coil_op={rho1,rho2,...,rhoN}
%
% where PhN have the same dimension as the sample voxel grid and rhoN are 
% spin states obtained from state() function.
%
%  parameters.*       - additional subfields may be required by your
%                       pulse sequence - check its documentation page
%
% Outputs:
%
%   This function returns whatever the pulse sequence returns.
%
% a.acharya@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=meshflow.m>

function answer=meshflow(spin_system,pulse_sequence,parameters)

% Check consistency
grumble(spin_system,pulse_sequence,parameters);

% Get diffusion and flow generator
F=flow_gen(spin_system,parameters);

% Get problem dimensions
spc_dim=spin_system.mesh.vor.ncells;  
spn_dim=size(spin_system.bas.basis,1); problem_dim=spc_dim*spn_dim;
report(spin_system,['lab space problem dimension     ' num2str(spc_dim)]);
report(spin_system,['spin space problem dimension    ' num2str(spn_dim)]);
report(spin_system,['Fokker-Planck problem dimension ' num2str(problem_dim)]);
parameters.spc_dim=spc_dim; parameters.spn_dim=spn_dim;

% Build Hamiltonian superoperator out of phantoms
polyadic_cores={{sparse([],[],[],spc_dim,spc_dim),...
                 sparse([],[],[],spn_dim,spn_dim)}};
for n=1:numel(parameters.H_ph)
    polyadic_cores{n}=cell(1,2);
    polyadic_cores{n}{1}=spdiags(parameters.H_ph{n}(:),0,spc_dim,spc_dim);
    polyadic_cores{n}{2}=parameters.H_op{n};
end
H=polyadic(polyadic_cores);

% Build relaxation superoperator out of phantoms
polyadic_cores={{sparse([],[],[],spc_dim,spc_dim),...
                 sparse([],[],[],spn_dim,spn_dim)}};
for n=1:numel(parameters.R_ph)
    polyadic_cores{n}=cell(1,2);
    polyadic_cores{n}{1}=spdiags(parameters.R_ph{n}(:),0,spc_dim,spc_dim);
    polyadic_cores{n}{2}=parameters.R_op{n};
end
R=polyadic(polyadic_cores);

% Build kinetics superoperator out of phantoms
polyadic_cores={{sparse([],[],[],spc_dim,spc_dim),...
                 sparse([],[],[],spn_dim,spn_dim)}};
for n=1:numel(parameters.K_ph)
    polyadic_cores{n}=cell(1,2);
    polyadic_cores{n}{1}=spdiags(parameters.K_ph{n}(:),0,spc_dim,spc_dim);
    polyadic_cores{n}{2}=parameters.K_op{n};
end
K=polyadic(polyadic_cores);

% Spin-independent spatial motion
F=polyadic({{F,opium(spn_dim,1)}});

% Dummy gradients for now
G=polyadic({{0}});

% Inflate polyadic objects if necessary
if ~ismember('polyadic',spin_system.sys.enable)
    H=inflate(H); R=inflate(R); K=inflate(K);
    G=inflate(G); F=inflate(F);
end

% Build the initial condition out of phantoms
parameters.rho0=0;
for n=1:numel(parameters.rho0_ph)
    parameters.rho0=parameters.rho0+kron(parameters.rho0_ph{n}(:),...
                                         parameters.rho0_st{n});
end

% Build the detection state out of phantoms
parameters.coil=0;
for n=1:numel(parameters.coil_ph)
    parameters.coil=parameters.coil+kron(parameters.coil_ph{n}(:),...
                                         parameters.coil_st{n});
end

% Call the pulse sequence
answer=pulse_sequence(spin_system,parameters,H,R,K,G,F);

end

% Consistency enforcement
function grumble(spin_system,pulse_sequence,parameters)
if ~isfield(spin_system,'mesh')
    error('mesh information is missing from the spin_system structure.');
end
if ~isfield(spin_system.mesh,'idx')
    error('indexing information is missing from spin_system.mesh structure.');
end
if ~isfield(spin_system.mesh,'vor')
    error('Voronoi tessellation information is missing from spin_system.mesh structure.');
end
if ~isa(pulse_sequence,'function_handle')
    error('pulse_sequence argument must be a function handle.');
end
if (~isfield(parameters,'H_ph'))||(~isfield(parameters,'H_op'))
    error('Hamiltonian phantom must be provided in parameters.H_ph and .H_op');
end
if (~isfield(parameters,'R_ph'))||(~isfield(parameters,'R_op'))
    error('relaxation phantom must be provided in parameters.R_ph and .R_op');
end
if (~isfield(parameters,'K_ph'))||(~isfield(parameters,'K_op'))
    error('kinetics phantom must be provided in parameters.K_ph and .K_op');
end
if (~isfield(parameters,'rho0_ph'))||(~isfield(parameters,'rho0_st'))
    error('initial state phantom must be provided in parameters.rho0_ph and .rho0_st');
end
if (~isfield(parameters,'coil_ph'))||(~isfield(parameters,'coil_st'))
    error('detection state phantom must be provided in parameters.coil_ph and .coil_st');
end
end

% Штык-нож на СВД имеет важное педагогическое значение. Если солдат или
% курсант на стрельбище ни разу не попал в мишень (или просто не выполнил
% норматив), он должен был примкнуть штык, ползком по-пластунски добрать-
% ся до мишени, заколоть ее штыком прямо в десятку и ползком же вернуться
% назад. После пары итераций меткость резко возрастала.
%
% Russian Internet folklore

