% Fokker-Planck imaging simulation context. Generates the Hamiltonian,
% the relaxation superoperator, the kinetics superoperator, the Fokker-
% Planck spatial dynamics generator (including diffusion and flow), gra-
% dient operators, and passes all of that to the pulse sequence, which
% should be supplied as a handle. Syntax:
%
%         answer=imaging(spin_system,pulse_sequence,parameters)
%
% Parameters:
%
%  pulse_sequence       -  pulse sequence function handle. See the
%                          experiments directory for the list of
%                          pulse sequences that ship with Spinach.
%
%   parameters.u       - X components of the velocity vectors
%                        for each point in the sample, m/s
%
%   parameters.v       - Y components of the velocity vectors
%                        for each point in the sample, m/s
%
%   parameters.w       - Z components of the velocity vectors
%                        for each point in the sample, m/s
%
%   parameters.diff    - diffusion coefficient or 3x3 tensor, m^2/s
%                        for situations when this parameter is the 
%                        same in every voxel
%
%   parameters.dxx     - Cartesian components of the diffusion
%   parameters.dxy       tensor for each voxel of the sample
%        ...
%   parameters.dzz
%
%   parameters.dims    - dimensions of the 3D box, meters
%
%   parameters.npts    - number of points in each dimension
%                        of the 3D box
%
%   parameters.deriv   - {'fourier'} uses Fourier diffe-
%                        rentiation matrices; {'period',n}
%                        requests n-point central finite-
%                        difference matrices with periodic
%                        boundary conditions
%
% Three types of phantoms must be specified. The relaxation theory phantom
% contains relaxation superoperators and their coefficients in each voxel,
% specified in the following way:
%
%                 parameters.rlx_ph={Ph1,Ph2,...,PhN}
%                 parameters.rlx_op={R1,R2,...,RN}
%
% where PhN have the same dimension as the sample voxel grid and RN are re-
% laxation superoperators. The initial condition phantom reflects the fact 
% that different voxels might start off in a different spin state. It must
% be specified in the following way:
%
%                 parameters.rho0_ph={Ph1,Ph2,...,PhN}
%                 parameters.rho0_op={rho1,rho2,...,rhoN}
%
% where PhN have the same dimension as the sample voxel grid and rhoN are 
% spin states obtained from state() function. The detection state phantom
% reflects the fact that different voxels might be detected at different 
% angles and with different sensitivity. It must be specified in the follo-
% wing way:
%
%                 parameters.coil_ph={Ph1,Ph2,...,PhN}
%                 parameters.coil_op={rho1,rho2,...,rhoN}
%
% where PhN have the same dimension as the sample voxel grid and rhoN are 
% spin states obtained from state() function.
%
% Outputs:
%
%  This function returns whatever the pulse sequence returns.
%
% Note: the direct product order is Z(x)Y(x)X(x)Spin, this cor-
%       responds to a column-wise vectorization of a 3D array
%       with dimensions ordered as [X Y Z].
%
% a.j.allami@soton.ac.uk
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=Imaging.m>

function answer=imaging(spin_system,pulse_sequence,parameters)

% Show the banner
banner(spin_system,'sequence_banner'); 

% Check consistency
grumble(spin_system,parameters);

% Set NMR assumptions
spin_system=assume(spin_system,'nmr');

% Call Spinach to build Hamiltonian superoperator
H=hamiltonian(spin_system);

% Process channel offsets
H=frqoffset(spin_system,H,parameters);

% Call Spinach to build kinetics superoperator
K=kinetics(spin_system);

% Get problem dimensions
spc_dim=prod(parameters.npts); spn_dim=size(H,1); problem_dim=spc_dim*spn_dim;
report(spin_system,['lab space problem dimension     ' num2str(spc_dim)]);
report(spin_system,['spin space problem dimension    ' num2str(spn_dim)]);
report(spin_system,['Fokker-Planck problem dimension ' num2str(problem_dim)]);
parameters.spc_dim=spc_dim; parameters.spn_dim=spn_dim;

% Build relaxation superoperator out of phantoms
polyadic_cores={{sparse([],[],[],spc_dim,spc_dim),...
                 sparse([],[],[],spn_dim,spn_dim)}};
for n=1:numel(parameters.rlx_ph)
    polyadic_cores{n}=cell(1,2);
    polyadic_cores{n}{1}=spdiags(parameters.rlx_ph{n}(:),0,spc_dim,spc_dim);
    polyadic_cores{n}{2}=parameters.rlx_op{n};
end
R=polyadic(polyadic_cores);

% Decide the initial state
if ~isfield(parameters,'rho0')

    % Build out of phantoms
    parameters.rho0=0;
    for n=1:numel(parameters.rho0_ph)
        parameters.rho0=parameters.rho0+kron(parameters.rho0_ph{n}(:),...
                                             parameters.rho0_st{n});
    end
    report(spin_system,'initial state built out of phantoms.');

else

    % Use initial state as received
    report(spin_system,'initial state supplied by the user.');

end

% Decide the coil state
if ~isfield(parameters,'coil')

    % Build out of phantoms
    parameters.coil=0;
    for n=1:numel(parameters.coil_ph)
        parameters.coil=parameters.coil+kron(parameters.coil_ph{n}(:),...
                                             parameters.coil_st{n});
    end
    report(spin_system,'coil state built out of phantoms.');

else

    % Use coil state as received
    report(spin_system,'coil state supplied by the user.');

end

% Put the same spin dynamics and chemistry into every voxel
H=polyadic({{opium(spc_dim,1),H}}); 
K=polyadic({{opium(spc_dim,1),K}});

% Get gradient operators
G=g2fplanck(spin_system,parameters);

% Get diffusion and hydrodynamics generators
F=v2fplanck(spin_system,parameters);

% Inflate polyadic objects
if ~ismember('polyadic',spin_system.sys.enable)
    H=complex(inflate(H)); R=complex(inflate(R)); 
    K=complex(inflate(K)); G=complex(inflate(G)); 
    F=complex(inflate(F));
end

% Call the pulse sequence
answer=pulse_sequence(spin_system,parameters,H,R,K,G,F);

end

% Consistency enforcement
function grumble(spin_system,parameters)

% Enforce relaxation phantom specification
if (~isfield(parameters,'rlx_ph'))||(~isfield(parameters,'rlx_op'))
    error('relaxation phantom must be provided in parameters.rlx_ph and parameters.rlx_op');
end
if (~iscell(parameters.rlx_ph))||(~iscell(parameters.rlx_op))
    error('parameters.rlx_ph and parameters.rlx_op must be cell arrays.');
end
if numel(parameters.rlx_ph)~=numel(parameters.rlx_op)
    error('parameters.rlx_ph and parameters.rlx_op must have the same number of elements.');
end

% Second dimension has one element if the sample is one-dimensional
if isscalar(parameters.npts), parameters.npts=[parameters.npts 1]; end

% If initial state is not provided, enforce building blocks
if ~isfield(parameters,'rho0')
    if (~isfield(parameters,'rho0_ph'))||(~isfield(parameters,'rho0_st'))
        error('initial state phantom must be provided in parameters.rho0_ph and parameters.rho0_st');
    end
    if (~iscell(parameters.rho0_ph))||(~iscell(parameters.rho0_st))
        error('parameters.rho0_ph and parameters.rho0_st must be cell arrays.');
    end
    if numel(parameters.rho0_ph)~=numel(parameters.rho0_st)
        error('parameters.rho0_ph and parameters.rho_st must have the same number of elements.');
    end
    for n=1:numel(parameters.rho0_ph)
        if ~all(size(parameters.rho0_ph{n})==parameters.npts)
            error('the size of all initial state phantoms must match parameters.npts');
        end
    end
    for n=1:numel(parameters.rlx_ph)
        for k=1:numel(parameters.rho0_ph)
            if ndims(parameters.rlx_ph{n})~=ndims(parameters.rho0_ph{k})
                error('all phantoms must have the same number of dimensions.');
            end
            if ~all(size(parameters.rlx_ph{n})==size(parameters.rho0_ph{k}))
                error('all phantoms must have the same size.');
            end
        end
    end
end

% If coil state is not provided, enforce building blocks
if ~isfield(parameters,'coil')
    if (~isfield(parameters,'coil_ph'))||(~isfield(parameters,'coil_st'))
        error('detection state phantom must be provided in parameters.coil_ph and parameters.coil_st');
    end
    if (~iscell(parameters.coil_ph))||(~iscell(parameters.coil_st))
        error('parameters.coil_ph and parameters.coil_st must be cell arrays.');
    end
    if numel(parameters.coil_ph)~=numel(parameters.coil_st)
        error('parameters.coil_ph and parameters.coil_st must have the same number of elements.');
    end
    for n=1:numel(parameters.coil_ph)
        if ~all(size(parameters.coil_ph{n})==parameters.npts)
            error('the size of all coil phantoms must match parameters.npts');
        end
    end
    for n=1:numel(parameters.rlx_ph)
        for m=1:numel(parameters.coil_ph)
            if ndims(parameters.rlx_ph{n})~=ndims(parameters.coil_ph{m})
                error('all phantoms must have the same number of dimensions.');
            end
            if ~all(size(parameters.rlx_ph{n})==size(parameters.coil_ph{m}))
                error('all phantoms must have the same size.');
            end
        end
    end
end

% Enforce phantom digitisation consistency
for n=1:numel(parameters.rlx_ph)
    if ~all(size(parameters.rlx_ph{n})==parameters.npts)
        error('the size of all relaxation phantoms must match parameters.npts');
    end
end

% Make sure pixel counts are integers
if isfield(parameters,'image_size')&&(~all(mod(parameters.image_size,2)))
    error('all elements of parameters.image_size must be odd integers.');
end

% Enforce no irrep mathematics
if ~isempty(spin_system.comp.sym_group)
    error('symmetry treatment is not supported in imaging simulations.');
end

end

% You have my sympathies with rotations theory - this is one of the 
% most treacherous parts of spin dynamics. At the recent ENC, a guy
% approached Malcolm Levitt and declared that all rotation conventi-
% ons in his book are off by a sign. I could see Malcolm's face tur-
% ning white. Three hours later, when I completed my poster reading
% cycle and returned to Malcolm's poster, the two were still there,
% arguing heatedly over a laptop with Mathematica running. Somebody
% asked me what they were discussing. I said, "religion".
%
% IK's email to Fred Mentink-Vigier, 2014

