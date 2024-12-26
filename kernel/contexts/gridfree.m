% Fokker-Planck magic angle spinning and SLE context. Generates a Liouvil-
% lian superoperator and passes it on to the pulse sequence function, which
% should be supplied as a handle. Syntax:
%
%     answer=gridfree(spin_system,pulse_sequence,parameters,assumptions)
%
% where pulse sequence is a function handle to one of the pulse sequences
% located in the experiments directory, assumptions is a string that would
% be passed to assume.m when the Hamiltonian is built and parameters is a
% structure with the following subfields:
%
%   parameters.rate     - spinning rate in Hz. Positive numbers
%                         for JEOL, negative for Varian and Bruker
%                         due to different rotation directions.
%
%   parameters.axis     - spinning axis, given as a normalized
%                         3-element vector
%
%   parameters.spins    - a cell array giving the spins that 
%                         the pulse sequence involves, e.g. 
%                         {'1H','13C'}
%
%   parameters.offset   - a cell array giving transmitter off-
%                         sets in Hz on each of the spins listed
%                         in parameters.spins array
%
%   parameters.max_rank - maximum D-function rank to retain in
%                         the solution (increase till conver-
%                         gence is achieved, approximately
%                         equal to the number of spinning si-
%                         debands in the spectrum)
%
%   parameters.tau_c - correlation times (in seconds) for rotational 
%                      diffusion. Single number for isotropic rotati-
%                      onal diffusion, and a 3x3 matrix for anisotro-
%                      pic rotational diffusion.
%
%   parameters.*       - additional subfields may be required by your
%                        pulse sequence - check its documentation page 
%
% The parameters structure is passed to the pulse sequence with the follo-
% wing additional parameters set:
%
%   parameters.spc_dim  - matrix dimension for the spatial
%                         dynamics subspace
%
%   parameters.spn_dim  - matrix dimension for the spin 
%                         dynamics subspace
%
% This function returns the powder average of whatever it is that the pulse
% sequence returns.
%
% Note: the choice of the Wigner D function rank truncation level depends on
%       on the spinning rate (the slower the spinning, the greater ranks are
%       required). 
%
% Note: rotational correlation times for SLE go into parameters.tau_c, not
%       inter.tau_c (the latter is only used by the Redfield theory module).
%
% Note: the state projector assumes a powder -- single crystal MAS is not
%       currently supported.
%
% Note: perturbative corrections to the rotating frame transformation are 
%       not supported - use singlerot.m instead.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=gridfree.m>

function answer=gridfree(spin_system,pulse_sequence,parameters,assumptions)

% Show the banner
banner(spin_system,'sequence_banner'); 

% Set common defaults
parameters=defaults(spin_system,parameters);

% Check consistency
grumble(spin_system,pulse_sequence,parameters,assumptions);

% Report to the user
report(spin_system,'building the Liouvillian...');

% Get spatial operators
[Lx,Ly,Lz,D,~]=sle_operators(parameters.max_rank);

% Set the assumptions
spin_system=assume(spin_system,assumptions);

% Get the Hamiltonian
[H,Q]=hamiltonian(spin_system);

% Disallow giant spin
if numel(Q)>2, error('giant spin model is not supported in this module.'); end

% Apply offsets
H=frqoffset(spin_system,H,parameters);

% Get relaxation and kinetics
R=relaxation(spin_system); 
K=kinetics(spin_system);

% Get problem dimensions
spc_dim=size(Lx,1); parameters.spc_dim=spc_dim;
spn_dim=size(H,1);  parameters.spn_dim=spn_dim;

% Inform the user
report(spin_system,['spin space problem dimension    ' num2str(spn_dim)]);
report(spin_system,['lab space space dimension       ' num2str(spc_dim)]);
report(spin_system,['Fokker-Planck problem dimension ' num2str(spn_dim*spc_dim)]);

% Project isotropic parts
H=polyadic({{speye(spc_dim),H}});
R=polyadic({{speye(spc_dim),R}});
K=polyadic({{speye(spc_dim),K}});

% Add anisotropic parts
for k=1:5
    for m=1:5
        H=H+polyadic({{D{k,m},Q{2}{k,m}}});
    end
end

% Rotation operator
if isfield(parameters,'rate')&&abs(parameters.rate)>0

    % Normalize axis vector
    spinning_axis=parameters.axis/norm(parameters.axis,2);
    
    % Report to the user
    report(spin_system,['spinning axis ort: ' num2str(spinning_axis)]);
    report(spin_system,['spinning rate: ' num2str(parameters.rate) ' Hz']);
    
    % Build rotation operator
    Hr=2*pi*parameters.rate*(spinning_axis(1)*Lx+spinning_axis(2)*Ly+spinning_axis(3)*Lz);
   
    % Add rotation operator
    H=H+polyadic({{Hr,speye(spn_dim)}});

end

% Diffusion operator
if isfield(parameters,'tau_c')
    
    % Build diffusion operator
    if isscalar(parameters.tau_c)
        
        % Isotropic rotational diffusion
        diff_iso=1/(6*parameters.tau_c);
        Rd=diff_iso*(Lx*Lx+Ly*Ly+Lz*Lz);
        
    elseif (size(parameters.tau_c,1)==3)&&...
           (size(parameters.tau_c,2)==3)
        
        % Anisotropic rotational diffusion
        dten=1./(6*parameters.tau_c);
        Rd=dten(1,1)*Lx*Lx+dten(1,2)*Lx*Ly+dten(1,3)*Lx*Lz+...
           dten(2,1)*Ly*Lx+dten(2,2)*Ly*Ly+dten(2,3)*Ly*Lz+...
           dten(3,1)*Lz*Lx+dten(3,2)*Lz*Ly+dten(3,3)*Lz*Lz;
        
    else
        
        % Complain and bomb out
        error('invalid correlation time specification.');
        
    end
    
    % Add diffusion operator
    R=R-polyadic({{Rd,speye(spn_dim)}});
    
end

% Inflate polyadics if necessary
if ~ismember('polyadic',spin_system.sys.enable)
    H=clean_up(spin_system,inflate(H),spin_system.tols.liouv_zero);
    R=clean_up(spin_system,inflate(R),spin_system.tols.liouv_zero);
    K=clean_up(spin_system,inflate(K),spin_system.tols.liouv_zero);
end

% Inform the user
H_whos=whos('H'); report(spin_system,['memory footprint of H array: ' num2str(H_whos.bytes/1024^3) ' GB']);
R_whos=whos('R'); report(spin_system,['memory footprint of R array: ' num2str(R_whos.bytes/1024^3) ' GB']);
K_whos=whos('K'); report(spin_system,['memory footprint of K array: ' num2str(K_whos.bytes/1024^3) ' GB']);

% Project initial and detection states into D[0,0,0]
P=spalloc(spc_dim,1,1); P(1)=1;
if isfield(parameters,'rho0')
    parameters.rho0=kron(P,parameters.rho0);
end
if isfield(parameters,'coil')
    parameters.coil=kron(P,parameters.coil);
end

% Report to the user
report(spin_system,'running the pulse sequence...');

% Call the pulse sequence
answer=pulse_sequence(spin_system,parameters,H,R,K);
    
end

% Default parameters
function parameters=defaults(spin_system,parameters)
if ~isfield(parameters,'offset')
    report(spin_system,'parameters.offset field not set, assuming zero offsets.');
    parameters.offset=zeros(size(parameters.spins));
end
if ~isfield(parameters,'verbose')
    report(spin_system,'parameters.verbose field not set, silencing array operations.');
    parameters.verbose=0;
end
end

% Consistency enforcement
function grumble(spin_system,pulse_sequence,parameters,assumptions)

% Rotating frames
if isfield(parameters,'rframes')
    error('numerical rotating frame transformation is not supported by SLE formalism.');
end

% Spherical grid
if isfield(parameters,'grid')
    error('spherical integration is a part of the SLE formalism, this parameter is unnecessary.');
end

% Formalism 
if ~ismember(spin_system.bas.formalism,{'zeeman-liouv','sphten-liouv'})
    error('this function is only available in Liouville space.');
end 

% Rotor rank
if ~isfield(parameters,'max_rank')
    error('parameters.max_rank subfield must be present.');
elseif (~isnumeric(parameters.max_rank))||(~isreal(parameters.max_rank))||...
       (mod(parameters.max_rank,1)~=0)||(parameters.max_rank<0)
    error('parameters.max_rank must be a positive real integer.');
end

% Spinning rate
if isfield(parameters,'rate')&&((~isnumeric(parameters.rate))||(~isreal(parameters.rate)))
    error('parameters.rate must be a real number.');
end

% Spinning axis
if isfield(parameters,'axis')&&((~isnumeric(parameters.axis))||(~isreal(parameters.axis))||...
                                (~isrow(parameters.axis))||(numel(parameters.axis)~=3))
    error('parameters.axis must be a row vector of three real numbers.');
end

% Check rotational correlation time
if isfield(parameters,'tau_c')
    if (~isreal(parameters.tau_c))||(~isnumeric(parameters.tau_c))
        error('parameters.tau_c must have real elements.');
    end
    if (numel(parameters.tau_c)~=1)&&(numel(parameters.tau_c)~=9)
        error('parameters.tau_c must be a scalar or a 3x3 matrix.');
    end
end

% Pulse sequence
if ~isa(pulse_sequence,'function_handle')
    error('pulse_sequence argument must be a function handle.');
end

% Assumptions
if ~ischar(assumptions)
    error('assumptions argument must be a character string.');
end

% Active spins
if isempty(parameters.spins)
    error('parameters.spins variable cannot be empty.');
elseif ~iscell(parameters.spins)
    error('parameters.spins variable must be a cell array.');
elseif ~all(cellfun(@ischar,parameters.spins))
    error('all elements of parameters.spins cell array must be strings.');
elseif any(~ismember(parameters.spins,spin_system.comp.isotopes))
    error('parameters.spins refers to a spin that is not present in the system.');
end

% Offsets
if isempty(parameters.offset)
    error('parameters.offset variable cannot be empty.');
elseif ~isnumeric(parameters.offset)
    error('parameters.offset variable must be an array of real numbers.');
elseif ~isfield(parameters,'spins')
    error('parameters.spins variable must be specified alongside parameters.offset variable.');
elseif numel(parameters.offset)~=numel(parameters.spins)
    error('parameters.offset variable must have the same number of elements as parameters.spins.');
end

end

% No god and no religion can survive ridicule. No political church, 
% no nobility, no royalty or other fraud, can face ridicule in a fair
% field, and live.
%
% Mark Twain

