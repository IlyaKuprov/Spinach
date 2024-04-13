% Hole burning experiment - a soft pulse follwed by a hard pi/2 
% observation pulse. The soft pulse is simulated using Fokker-
% Planck formalism. Syntax:
%
%          fid=holeburn(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     parameters.pulse_frq - frequency of the soft pulse, Hz
%
%     parameters.pulse_phi - phase of the soft pulse, rad
%
%     parameters.pulse_pwr - power of the soft pulse, rad/s
%
%     parameters.pulse_dur - duration of the sof pulse, s
%
%     parameters.pulse_rnk - Fokker-Planck cut-off rank
%
%     parameters.offset    - receiver offset for the time
%                            domain detection, Hz
%
%     parameters.sweep     - sweep width for time domain
%                            detection, Hz
%
%     parameters.npoints   - number of points in the free
%                            induction decay 
%
%     parameters.rho0      - initial state
%
%     parameters.coil      - detection state
%
%     parameters.method    - soft puse propagation method,
%                            'expv' for Krylov propagation,
%                            'expm' for exponential propa-
%                            gation, 'evolution' for Spin-
%                            ach evolution function
%
%    H     - Hamiltonian matrix, received from context 
%            function
%
%    R     - relaxation superoperator, received from 
%            context function
%
%    K     - kinetics superoperator, received from 
%            context function
%
% Output:
%
%     fid  -  free induction decay seen by the state speci-
%             fied in parameters coil after the hole burning 
%             pulse followed by a hard pi/2 pulse. 
%
% Note: the rank parameter should be increased until conver-
%       gence is achieved in the output.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=holeburn.m>

function fid=holeburn(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Pulse operators
Ep=operator(spin_system,'L+',parameters.spins{1});
Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i;

% A soft pulse
parameters.pulse_frq=parameters.pulse_frq-parameters.offset;
rho=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0,parameters.pulse_frq,...
                                        parameters.pulse_pwr,parameters.pulse_dur,...
                                        parameters.pulse_phi,parameters.pulse_rnk,...
                                        parameters.method);

% A hard pulse 
parameters.rho0=step(spin_system,Ey,rho,pi/2);
                                  
% Acquisition
fid=acquire(spin_system,parameters,H,R,K);

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function is only available in Liouville space.');
end
if ~isfield(parameters,'pulse_frq')
    error('pulse frequency must be specified in parameters.pulse_frq field.');
end
if (~isnumeric(parameters.pulse_frq))||(~isreal(parameters.pulse_frq))
    error('parameters.pulse_frq must be an array of reals.');
end
if ~isfield(parameters,'pulse_pwr')
    error('pulse power must be specified in parameters.pulse_pwr field.');
end
if (~isnumeric(parameters.pulse_pwr))||(~isreal(parameters.pulse_pwr))
    error('parameters.pulse_pwr must be an array of reals.');
end
if ~isfield(parameters,'pulse_dur')
    error('pulse duration must be specified in parameters.pulse_dur field.');
end
if (~isnumeric(parameters.pulse_dur))||(~isreal(parameters.pulse_dur))
    error('parameters.pulse_dur must be an array of reals.');
end
if ~isfield(parameters,'pulse_phi')
    error('pulse phase must be specified in parameters.pulse_phi field.');
end
if (~isnumeric(parameters.pulse_phi))||(~isreal(parameters.pulse_phi))||...
   (numel(parameters.pulse_phi)~=1)
    error('parameters.pulse_phi must be a real scalar.');
end
if ~isfield(parameters,'pulse_rnk')
    error('pulse grid rank must be specified in parameters.pulse_rnk field.');
end
if (~isnumeric(parameters.pulse_rnk))||(~isreal(parameters.pulse_rnk))||...
   (numel(parameters.pulse_rnk)~=1)||(mod(parameters.pulse_rnk,1)~=0)
    error('parameters.pulse_rnk must be a real integer.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'method')
    error('shaped pulse simulation method must be specified in parameters.method field.');
end
if ~isfield(parameters,'sweep')
    error('width of the detection window must be specified in parameters.sweep field.');
end
if (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
   (~isscalar(parameters.sweep))||(parameters.sweep<=0)
    error('parameters.sweep must be a positive real scalar.');
end
if ~isfield(parameters,'npoints')
    error('number of points in the FID must be specified in parameters.npoints field.');
end
if (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
   (~isscalar(parameters.npoints))||(parameters.npoints<1)||...
   (mod(parameters.npoints,1)~=0)
    error('parameters.npoints must be a positive real integer.');
end
end

% Much of the social history of the Western world, over the past
% three decades, has been a history of replacing what worked with
% what sounded good.
%
% Thomas Sowell

