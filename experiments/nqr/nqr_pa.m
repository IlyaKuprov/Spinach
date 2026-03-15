% Nuclear quadrupole resonance soft pulse-acquire experiment. Idealised ac-
% quisition with infinite bandwidth is done. Syntax:
%
%              spectrum=nqr_pa(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep           vector with two elements giving
%                               the spectrum window extents, Hz
%
%    parameters.npoints         number of points in the spectrum
%
%    parameters.rho0            initial state
%
%    parameters.coil            detection state
%
%    parameters.Lx              Lx and Ly operators that go into
%    parameters.Ly              the RF Hamiltonian
%                               
%    parameters.rf_frq          RF irradiation frequency, Hz
%
%    parameters.rf_pwr          the multiplier (rad/s) in front 
%                               of [Lx*cos(ωt)+Ly*sin(ωt)] in the
%                               RF Hamiltonian
%
%    parameters.rf_dur          pulse duration, seconds
%
%    H  - Hamiltonian matrix, received from context function
%
%    R  - relaxation superoperator, received from context function
%
%    K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    spectrum  - the spectrum of the system with the specified
%                starting state detected on the specified coil
%                state within the frequency interval requested
%
% Note: relaxation must be present in the system dynamics, or the 
%       matrix inversion operation would fail to converge. The re-
%       laxation matrix R should *not* be thermalised.
%
% ilya.kuprov@weizmann.ac.il
% lewis.robertson@csiro.au
%
% <https://spindynamics.org/wiki/index.php?title=nqr_pa.m>

function spectrum=nqr_pa(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Project pulse operators
Lx=kron(speye(parameters.spc_dim),parameters.Lx);
Ly=kron(speye(parameters.spc_dim),parameters.Ly);

% Apply the soft off-resonance pulse
parameters.rho0=shaped_pulse_af(spin_system,H+1i*R+1i*K,Lx,Ly,parameters.rho0,...
                                parameters.rf_frq,parameters.rf_pwr,...
                                parameters.rf_dur,0,2,'expm');
        
% Call frequency domain acquisition
spectrum=slowpass(spin_system,parameters,H,R,K);

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
if ~isfield(parameters,'sweep')
    error('sweep window must be specified in parameters.sweep variable.');
end
if (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
   (~isrow(parameters.sweep))||(numel(parameters.sweep)~=2)||...
   (parameters.sweep(1)>=parameters.sweep(2))
    error('parameters.sweep must be a row vector with two elements in ascending order.');
end
if ~isfield(parameters,'npoints')
    error('point count must be specified in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(numel(parameters.npoints)~=1)||...
   (~isreal(parameters.npoints))||(parameters.npoints<1)||(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must be a positive integer.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'Lx')
    error('RF Hamiltonian operator Lx must be specified in parameters.Lx variable.');
end
if (~isnumeric(parameters.Lx))||(~ismatrix(parameters.Lx))
    error('parameters.Lx must be a numeric matrix.');
end
if ~isfield(parameters,'Ly')
    error('RF Hamiltonian operator Ly must be specified in parameters.Ly variable.');
end
if (~isnumeric(parameters.Ly))||(~ismatrix(parameters.Ly))
    error('parameters.Ly must be a numeric matrix.');
end
if ~isequal(size(parameters.Lx),size(parameters.Ly))
    error('parameters.Lx and parameters.Ly must have the same dimension.');
end
if ~isfield(parameters,'rf_frq')
    error('RF irradiation frequency must be specified in parameters.rf_frq variable.');
end
if (~isnumeric(parameters.rf_frq))||(~isreal(parameters.rf_frq))||(~isscalar(parameters.rf_frq))
    error('parameters.rf_frq must be a real scalar.');
end
if ~isfield(parameters,'rf_pwr')
    error('RF irradiation power must be specified in parameters.rf_pwr variable.');
end
if (~isnumeric(parameters.rf_pwr))||(~isreal(parameters.rf_pwr))||...
   (~isscalar(parameters.rf_pwr))
    error('parameters.rf_pwr must be a real scalar.');
end
if ~isfield(parameters,'rf_dur')
    error('RF pulse duration must be specified in parameters.rf_dur variable.');
end
if (~isnumeric(parameters.rf_dur))||(~isreal(parameters.rf_dur))||...
   (~isscalar(parameters.rf_dur))||(parameters.rf_dur<0)
    error('parameters.rf_dur must be a non-negative real scalar.');
end
end

% Except where the backbone is an actual backbone. Ever been 
% to Malacath's realm...? Nasty stuff. But, back to the busi-
% ness at hand.
%
% Sheogorath

