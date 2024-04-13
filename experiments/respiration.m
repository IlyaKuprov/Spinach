% RESPIRATION cross-polarisation method described in the paper from
% the Aarhus group (http://dx.doi.org/10.1021/jz3000905). Syntax:
%
%           fid=respiration(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep              sweep width, Hz
%
%    parameters.npoints            number of points in the FID
%
%    parameters.rho0               initial state
%
%    parameters.coil               detection state
%
%    parameters.nloops             number of RESPIRATION loops
%
%    parameters.theta              the angle of the ideal pulse
%                                  at the end of each loop
%
%    parameters.spins              working spins, e.g. {'1H',13C'}
%
%    H     - Hamiltonian matrix, received from context function
%
%    R     - relaxation superoperator, received from context function
%
%    K     - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid   - free induction decay as seen by the state specified
%            in parameters parameters.coil
%
% venkata-subbarao.redrouthu@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=respiration.m>

function fid=respiration(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Generate and project pulse operators
Hp=operator(spin_system,'L+',parameters.spins{1});
Cp=operator(spin_system,'L+',parameters.spins{2});
Hp=kron(speye(parameters.spc_dim),Hp);
Cp=kron(speye(parameters.spc_dim),Cp);
Hx=(Hp+Hp')/2; Cx=(Cp+Cp')/2; 

% Build the Liouvillian
L=H+1i*R+1i*K;

% Initial condition
rho=parameters.rho0;

% RESPIRATION loop
for n=1:parameters.nloops
    
    % Update the user
    report(spin_system,['RESPIRATION loop ' num2str(n) '/'...
                        num2str(parameters.nloops) '...']);
    
    % +X pulse on protons
    rho=step(spin_system,L+2*pi*2*parameters.rate*Hx,rho,1/(2*parameters.rate));
    
    % -X pulse on protons
    rho=step(spin_system,L-2*pi*2*parameters.rate*Hx,rho,1/(2*parameters.rate));
    
    % Ideal theta pulse on both
    rho=step(spin_system,Cx+Hx,rho,parameters.theta);
    
end

% Acquisition
[L,rho]=decouple(spin_system,L,rho,{'1H'});
fid=evolution(spin_system,L,parameters.coil,rho,...
              1/parameters.sweep,parameters.npoints-1,'observable');
          
end

% Consistency enforcement
function grumble(parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
end
if (~isnumeric(parameters.sweep))||(numel(parameters.sweep)~=1)||...
   (~isreal(parameters.sweep))||(parameters.sweep<=0)
    error('parameters.sweep should be a positive real number.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(numel(parameters.npoints)~=1)||...
   (~isreal(parameters.npoints))||(parameters.npoints<1)||...
   (mod(parameters.npoints,1)~=0)
    error('parameters.npoints should be a positive integer.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'nloops')
    error('the number of loops must be specified in parameters.nloops variable.');
end
if ~isfield(parameters,'theta')
    error('short pulse angle must be specified in parameters.theta variable.');
end
end

% If people do not believe that mathematics is simple, 
% it is only because they do not realise how complica-
% ted life is.
% 
% John von Neumann

