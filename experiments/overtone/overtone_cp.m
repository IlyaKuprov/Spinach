% Cross-polarization overtone experiment. Syntax:
%
%           spectrum=overtone_cp(spin_system,parameters,H,R,K)
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
%    parameters.Nx              X Zeeman operator on the
%                               quadrupolar nucleus
%
%    parameters.Hx              X Zeeman operator on the
%                               spin-1/2 nucleus
%
%    parameters.rf_frq          spin-lock frequency offset from
%                               the overtone frequency on the 
%                               quadrupolar nucleus, Hz
%
%    parameters.rf_pwr          a vector of spin-lock powers on
%                               the quadrupolar nucleus (first
%                               element) and the spin-1/2 nucleus
%                               (second element), rad/s
%
%    parameters.rf_dur          spin-lock pulse duration, seconds
%
% Relaxation must be present in the system dynamics, or the matrix
% inversion operation below would fail to converge. The relaxation
% matrix should *not* be thermalized.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=Overtone_cp.m>

function spectrum=overtone_cp(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Get the overtone frequency
ovt_frq=-2*spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi);

% Project pulse operators
Nx=kron(speye(parameters.spc_dim),parameters.Nx);
Hx=kron(speye(parameters.spc_dim),parameters.Hx);

% Choose the method
switch parameters.method
    
    case 'average'

        % Get the frequency
        omega=2*pi*ovt_frq-2*pi*parameters.rf_frq;
        
        % Get the pulse Hamiltonian
        pulseop=average(spin_system,parameters.rf_pwr(1)*Nx/2,H+parameters.rf_pwr(2)*Hx,...
                                    parameters.rf_pwr(1)*Nx/2,omega,'matrix_log');

        % Apply the pulse
        parameters.rho0=propagator(spin_system,pulseop,parameters.rf_dur)*parameters.rho0;
        
    case 'fplanck'
        
        % Get the frequency
        omega=ovt_frq-parameters.rf_frq;

        % Apply the pulse
        parameters.rho0=shaped_pulse_af(spin_system,H+parameters.rf_pwr(2)*Hx,Nx,0*Nx,...
                                        parameters.rho0,omega,parameters.rf_pwr(1),...
                                        parameters.rf_dur,-pi/2,2,'expm');
        
end

% Call the acquisition
spectrum=overtone_a(spin_system,parameters,H,R,K);

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
if ~isfield(parameters,'sweep')
    error('spectral range must be specified in parameters.sweep variable.');
end
if (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
   (numel(parameters.sweep)~=2)
    error('parameters.sweep vector must have two real elements.');
end
if ~isfield(parameters,'Nx')
    error('quadrupolar nucleus X pulse operator must be supplied in parameters.Nx field.');
end
if ~isfield(parameters,'Hx')
    error('spin-1/2 nucleus X pulse operator must be supplied in parameters.Hx field.');
end
if ~isfield(parameters,'rf_frq')
    error('parameters.rf_frq field is missing.');
end
if ~isfield(parameters,'rf_pwr')
    error('parameters.rf_pwr field is missing.');
end
if ~isfield(parameters,'rf_dur')
    error('parameters.rf_dur field is missing.');
end
end

% Violence does not stem from a psychopathic lack of 
% morality. Quite the reverse: it comes from the ex-
% ercise of perceived moral rights and obligations.
% 
% Tage Rai

