% Overtone signal acquisition experiment in the frequency domain. Syntax:
%
%           spectrum=overtone_a(spin_system,parameters,H,R,K)
%
% where H is the hamiltonian commutation superoperator, R is the relaxation
% superoperator and K is the chemical kinetics superoperator. The following
% parameters are required:
%
%    parameters.sweep           vector with two elements giving
%                               the spectrum frequency extents
%                               in Hz around the overtone frequency
%
%    parameters.npoints         number of points in the spectrum
%
%    parameters.rho0            initial state
%
%    parameters.coil            detection state
%
% Relaxation must be present in the system dynamics, or the matrix
% inversion operation below would fail to converge. The relaxation
% matrix should *not* be thermalized.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=Overtone_a.m>

function spectrum=overtone_a(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Get the overtone frequency
ovt_frq=-2*spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi);

% Convert sweep to absolute frequencies
parameters.sweep=ovt_frq-parameters.sweep;

% Call slowpass
spectrum=slowpass(spin_system,parameters,H,R,K);

end

% Consistency enforcement
function grumble(parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(numel(parameters.npoints)~=1)||...
   (~isreal(parameters.npoints))||(parameters.npoints<1)||(mod(parameters.npoints,1)~=0)
    error('parameters.npoints should be a positive integer.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'sweep')
    error('spectral range must be specified in parameters.sweep variable.');
end
if (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||(numel(parameters.sweep)~=2)
    error('parameters.sweep vector must have two real elements');
end
end

% You would never earn a lot of money if you
% think that money is earned. 
%
% A Russian saying

