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
%    parameters.spins           overtone-active nucleus, specified as a
%                               single-element cell array
%
%    parameters.spc_dim         Fokker-Planck spatial dimension
%
%    parameters.rho0            initial state
%
%    parameters.coil            detection state
%
%    parameters.Lx              X Zeeman operator on the
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
% Outputs:
%
%    spectrum                   overtone spectrum
%
% Relaxation must be present in the system dynamics, or the matrix
% inversion operation below would fail to converge. The relaxation
% matrix should *not* be thermalized.
%
% ilya.kuprov@weizmann.ac.il
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
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'spins')
    error('overtone-active nucleus must be specified in parameters.spins variable.');
elseif (~iscell(parameters.spins))||(numel(parameters.spins)~=1)
    error('parameters.spins should be a single-element cell array.');
end
if ~isfield(parameters,'spc_dim')
    error('Fokker-Planck dimension should be specified in parameters.spc_dim variable.');
elseif (~isnumeric(parameters.spc_dim))||(~isreal(parameters.spc_dim))||...
       (numel(parameters.spc_dim)~=1)||(parameters.spc_dim<1)||...
       (mod(parameters.spc_dim,1)~=0)
    error('parameters.spc_dim should be a positive integer.');
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
elseif (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
       (numel(parameters.sweep)~=2)
    error('parameters.sweep vector must have two real elements.');
end
if ~isfield(parameters,'Lx')
    error('quadrupolar nucleus X pulse operator must be supplied in parameters.Lx field.');
end
if ~isfield(parameters,'rf_frq')
    error('parameters.rf_frq field is missing.');
elseif (~isnumeric(parameters.rf_frq))||(~isreal(parameters.rf_frq))||...
       (numel(parameters.rf_frq)~=1)
    error('parameters.rf_frq should be a real scalar.');
end
if ~isfield(parameters,'rf_pwr')
    error('parameters.rf_pwr field is missing.');
elseif (~isnumeric(parameters.rf_pwr))||(~isreal(parameters.rf_pwr))||...
       (numel(parameters.rf_pwr)~=1)
    error('parameters.rf_pwr should be a real scalar.');
end
if ~isfield(parameters,'rf_dur')
    error('parameters.rf_dur field is missing.');
elseif (~isnumeric(parameters.rf_dur))||(~isreal(parameters.rf_dur))||...
       (numel(parameters.rf_dur)~=1)||(parameters.rf_dur<=0)
    error('parameters.rf_dur should be a positive real scalar.');
end
if ~isfield(parameters,'method')
    error('pulse simulation method must be specified in parameters.method variable.');
elseif ~ismember(parameters.method,{'average','fplanck'})
    error('parameters.method must be ''average'' or ''fplanck''.');
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

