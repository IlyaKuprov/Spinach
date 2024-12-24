% Soft pulse followed by acquisition. The soft pulse is simulated 
% using the Fokker-Planck formalism. Syntax:
%
%          fid=sp_acquire(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     parameters.pulse_frq - frequency of the soft pulse,
%                            relative to the frequency of 
%                            the current rotating frame, Hz
%
%     parameters.pulse_phi - phase of the soft pulse, rad
%
%     parameters.pulse_pwr - power of the soft pulse, rad/s
%
%     parameters.pulse_dur - duration of the sof pulse, s
%
%     parameters.pulse_rnk - Fokker-Planck cut-off rank,
%                            a small integer: start with 2
%                            and increase until the answer
%                            stops changing 
%
%     parameters.offset    - transmitter / receiver offset
%                            for the time domain pulses and
%                            detection, relative to the cur-
%                            rent rotating frame, Hz
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
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=sp_acquire.m>

function fid=sp_acquire(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Pulse operators
Ep=operator(spin_system,'L+',parameters.spins{1});
Ep=kron(speye(parameters.spc_dim),Ep);
Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i;

% A soft pulse
parameters.pulse_frq=parameters.pulse_frq-parameters.offset;
parameters.rho0=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0,...
                                parameters.pulse_frq,parameters.pulse_pwr,...
                                parameters.pulse_dur,parameters.pulse_phi,...
                                parameters.pulse_rnk,parameters.method);
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
if (~isnumeric(parameters.pulse_frq))||(~isreal(parameters.pulse_frq))||...
   (numel(parameters.pulse_frq)~=1)
    error('parameters.pulse_frq must be a real scalar.');
end
if ~isfield(parameters,'pulse_pwr')
    error('pulse power must be specified in parameters.pulse_pwr field.');
end
if (~isnumeric(parameters.pulse_pwr))||(~isreal(parameters.pulse_pwr))||...
   (numel(parameters.pulse_pwr)~=1)||(parameters.pulse_pwr<0)
    error('parameters.pulse_pwr must be a non-negative real scalar.');
end
if ~isfield(parameters,'pulse_dur')
    error('pulse duration must be specified in parameters.pulse_dur field.');
end
if (~isnumeric(parameters.pulse_dur))||(~isreal(parameters.pulse_dur))||...
   (numel(parameters.pulse_dur)~=1)||(parameters.pulse_dur<=0)
    error('parameters.pulse_dur must be a positive real scalar.');
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

% It's now very common to hear people say "I'm rather offended by that"
% As if that gives them certain rights. It's actually nothing more [...]
% than a whine. "I find that offensive" -- it has no meaning; it has no
% purpose; it has no reason to be respected as a phrase. "I am offended
% by that" - well, so fucking what.
%
% Stephen Fry

