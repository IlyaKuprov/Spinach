% Phase-sensitive homonuclear NOESY pulse sequence. Syntax:
%
%             fid=noesy(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep     sweep widths, Hz
%
%    parameters.npoints   number of points for both dimensions
%
%    parameters.spins     nuclei on which the sequence runs,
%                         specified as '1H', '13C', etc.
%
%    parameters.tmix      mixing time, seconds
%
%    parameters.decouple  spins to be decoupled, specified either
%                         by name, e.g. {'13C','1H'}, or by a list
%                         of numbers, e.g. [1 2]
%
%    parameters.rho0      initial state; skip this and specify
%                         parameters.needs={'rho_eq'} to start
%                         from exact thermal equilibrium
%
%    H - Hamiltonian matrix, received from context function
%
%    R - relaxation superoperator, received from context function
%
%    K - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid.cos,fid.sin - two components of the FID for F1 hyper-
%                      complex processing
%
% Note: this function is used for extreme simulations (proteins
%       and nucleic acids) - its layout is optimised for minimum
%       memory footprint rather than CPU time.
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
% hannah.hogben@chem.ox.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=noesy.m>

function fid=noesy(spin_system,parameters,H,R,K)

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Coherent evolution timestep
timestep=1./parameters.sweep;

% Detection state
coil=state(spin_system,'L+',parameters.spins{1},'cheap');

% Pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx=(Lp+Lp')/2; Ly=(Lp-Lp')/2i;

% Decoupling
if isfield(parameters,'decouple')
    [H,parameters.rho0]=decouple(spin_system,H,parameters.rho0,...
                                               parameters.decouple);
end

% First pulse
rho=step(spin_system,Lx,parameters.rho0,pi/2);

% Phase cycle specification
Op2={Lx,Ly,Lx,Ly}; An2={+pi/2,+pi/2,-pi/2,-pi/2};
Op3={Ly,Ly,Ly,Ly}; An3={+pi/2,+pi/2,+pi/2,+pi/2};

% FID phase cycle
fids=cell(1,4);

% Phase cycle loop
for n=1:4

    % F1 evolution
    rho_stack=evolution(spin_system,H+1i*R+1i*K,[],rho,timestep(1),...
                        parameters.npoints(1)-1,'trajectory');
    % Second pulse
    rho_stack=step(spin_system,Op2{n},rho_stack,An2{n});

    % Homospoil
    rho_stack=homospoil(spin_system,rho_stack,'destroy');

    % Mixing time
    rho_stack=evolution(spin_system,1i*R+1i*K,[],...
                        rho_stack,parameters.tmix,1,'final');
    % Homospoil
    rho_stack=homospoil(spin_system,rho_stack,'destroy');

    % Third pulse
    rho_stack=step(spin_system,Op3{n},rho_stack,An3{n});
    
    % F2 evolution and detection
    fids{n}=evolution(spin_system,H+1i*R+1i*K,coil,rho_stack,...
                      timestep(2),parameters.npoints(2)-1,'observable');
                  
end

% Axial peak elimination
fid.cos=fids{1}-fids{3}; fid.sin=fids{2}-fids{4};

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function is only available for sphten-liouv and zeeman-liouv formalisms.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=2
    error('parameters.sweep array should have exactly two elements.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=1
    error('parameters.spins cell array should have exactly one element.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=2
    error('parameters.npoints array should have exactly two elements.');
end
if ~isfield(parameters,'tmix')
    error('mixing time should be specified in parameters.tmix variable.');
elseif numel(parameters.tmix)~=1
    error('parameters.tmix array should have exactly one element.');
end
end

% According to a trade legend, Anil Kumar had to run the very first
% NOESY experiment on a Saturday -- his supervisors viewed it as a
% waste of valuable spectrometer time.

