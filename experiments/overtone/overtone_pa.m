% Overtone soft pulse-acquire experiment. Syntax:
%
%           spectrum=overtone_pa(spin_system,parameters,H,R,K)
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
%    parameters.Lx
%    parameters.Ly              X and Y Zeeman operators on the
%                               quadrupolar nucleus
%
%    parameters.rf_frq          pulse frequency offset from
%                               the overtone frequency on the 
%                               quadrupolar nucleus, Hz
%
%    parameters.rf_pwr          pulse power on the quadrupolar 
%                               nucleus, rad/s
%
%    parameters.rf_dur          pulse duration, seconds
%
%    parameters.method          'average' uses the average Hamil-
%                               tonian theory, 'fplanck' uses
%                               Fokker-Planck formalism
%
% Relaxation must be present in the system dynamics, or the matrix
% inversion operation below would fail to converge. The relaxation
% matrix should *not* be thermalized.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=Overtone_pa.m>

function spectrum=overtone_pa(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Get the overtone frequency
ovt_frq=-2*spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi);

% Project pulse operators
Lx=kron(speye(parameters.spc_dim),parameters.Lx);

% Choose the method
switch parameters.method
    
    case 'average'

        % Get the pulse frequency
        omega=2*pi*ovt_frq-2*pi*parameters.rf_frq;
        
        % Get the pulse Hamiltonian
        pulseop=parameters.rf_pwr*Lx;
        pulseop=average(spin_system,pulseop/2,H,pulseop/2,omega,'matrix_log');

        % Apply the pulse
        parameters.rho0=propagator(spin_system,pulseop,parameters.rf_dur)*parameters.rho0;
        
    case 'fplanck'
        
        % Get the pulse frequency
        pulse_frq=ovt_frq-parameters.rf_frq;

        % Apply the pulse
        parameters.rho0=shaped_pulse_af(spin_system,H,Lx,0*Lx,parameters.rho0,pulse_frq,...
                                        parameters.rf_pwr,parameters.rf_dur,-pi/2,2,'expm');
        
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
   (~isreal(parameters.npoints))||(parameters.npoints<1)||(mod(parameters.npoints,1)~=0)
    error('parameters.npoints should be a positive integer.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
end

% The Nirvana Fallacy is a name given to the informal fallacy of
% comparing actual things with unrealistic, idealized alternati-
% ves. [...] By creating a false dichotomy that presents one op-
% tion which is obviously advantageous - while at the same time
% being completely implausible - a person using the nirvana fal-
% lacy can attack any opposing idea because it is imperfect. Un-
% der this fallacy, the choice is not between real world soluti-
% ons; it is, rather, a choice between one realistic achievable
% possibility and another unrealistic solution that could in so-
% me way be "better".
%
% Wikipedia

