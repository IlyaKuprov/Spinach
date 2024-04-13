% Budker group style gamma-weighted pulse-acquire sequence in zero 
% field. Uses gamma-weighted initial state (corresponding to using 
% a pre-polarisation magnet at high temperature), gamma-weighted
% pulse operators and gamma-weighted detection state. Syntax:
%
%            fid=zerofield(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep - the width of the spectral window (Hz)
%
%    parameters.npoints  - number time steps in the simulation
%
%    parameters.detection - 'uniaxial' to emulate common ZULF
%                           hardware, 'quadrature' for proper
%                           frequency sign discrimination
%
%    parameters.flip_angle - pulse flip angle in radians for
%                            protons; for other nuclei, this
%                            will be scaled by the gamma ratio
%
%    H  - Hamiltonian matrix, received from context function
%
%    R  - relaxation superoperator, received from context function
%
%    K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid - free induction decay
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=zerofield.m>

function fid=zerofield(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Magnetogyric ratio weights relative to 1H
weights=spin_system.inter.gammas/spin('1H');

% Get gamma-weighted initial state
rho=sparse(0);
for n=1:spin_system.comp.nspins
    rho=rho+weights(n)*state(spin_system,{'Lz'},{n});
end

% Get gamma-weighted detection state
coil=sparse(0);
for n=1:spin_system.comp.nspins
    coil=coil+weights(n)*state(spin_system,{'L+'},{n});
end

% Get gamma-weighted pulse operator
Sy=sparse(0);
for n=1:spin_system.comp.nspins
    Sy=Sy+weights(n)*operator(spin_system,{'Ly'},{n});
end

% Apply the pulse
rho=step(spin_system,Sy,rho,parameters.flip_angle);

% Compute the digitisation parameters
timestep=1/parameters.sweep;

% Run the simulation
fid=evolution(spin_system,L,coil,rho,timestep,...
              parameters.npoints-1,'observable');
          
% Emulate detection hardware
switch parameters.detection
    
    case 'quadrature'
        
        % Do nothing
        
    case 'uniaxial'
        
        % Destroy imaginary part
        fid=real(fid);
        
end

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
if ~isfield(parameters,'detection')
    error('detection mode must be specified in parameters.detection');
end
if ~ismember(parameters.detection,{'uniaxial','quadrature'})
    error('parameters.detection must be ''uniaxial'' or ''quadrature''');
end
if ~isfield(parameters,'flip_angle')
    error('proton filip angle must be specified in parameters.flip_angle');
end
end

% The only real failure in life is not to 
% be true to the best one knows.
%
% Buddha

