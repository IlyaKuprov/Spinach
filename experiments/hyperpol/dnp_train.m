% Repeated piecewise-constant microwave block for pulsed DNP. Syntax:
%
%           signal=dnp_train(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.rho0        - initial Liouville-space column vector
%
%    parameters.coil        - detection column vector
%
%    parameters.pulse_ops   - {Ex,Ey} control superoperators
%
%    parameters.amplitudes  - {x,y} Cartesian row vectors in rad/s
%
%    parameters.slice_durs  - row vector of positive durations (s)
%
%    parameters.nloops      - positive integer number of blocks
%
%    H, R, K                - Hamiltonian, relaxation, and kinetics
%                             superoperators supplied by the context
%
% Outputs:
%
%    signal - complex column with nloops+1 entries: coil overlap
%             before irradiation and after each complete block
%
% Use a static Liouville-space context such as powder or crystal.
% A single block propagator is reused; the drift and RF pattern
% must be the same on every repeat. No normalisation is applied.
%
% Talos <talos@spindynamics.org>

function signal=dnp_train(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Propagate the identity through one block without storing slice states
L=H+1i*R+1i*K;
P=shaped_pulse_xy(spin_system,L,parameters.pulse_ops,parameters.amplitudes,...
                  parameters.slice_durs,speye(size(L)),'expm-pwc');

% Detect at the start and after every repeat
rho=parameters.rho0;
signal=zeros(parameters.nloops+1,1,'like',1i);
signal(1)=parameters.coil'*rho;
for n=1:parameters.nloops

    % Advance one period and detect
    rho=P*rho;
    signal(n+1)=parameters.coil'*rho;

end

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('dnp_train requires a Liouville-space basis.');
end
dim=spin_system.bas.offsets(end);
if ~all(cellfun(@(M)isnumeric(M)&&isequal(size(M),[dim dim]),{H,R,K}))
    error('H, R, and K must be square matrices matching the spin basis.');
end
if ~isnumeric(parameters.rho0)||~isequal(size(parameters.rho0),[dim 1])||...
   ~isnumeric(parameters.coil)||~isequal(size(parameters.coil),[dim 1])
    error('parameters.rho0 and parameters.coil must be spin basis column vectors.');
end
if ~iscell(parameters.pulse_ops)||~isequal(size(parameters.pulse_ops),[1 2])||...
   ~all(cellfun(@(M)isnumeric(M)&&isequal(size(M),[dim dim]),parameters.pulse_ops))
    error('parameters.pulse_ops must be {Ex,Ey} with operators matching the spin basis.');
end
validateattributes(parameters.slice_durs,{'double','single'},{'row','nonempty','real','finite','positive'});
validateattributes(parameters.nloops,{'numeric'},{'scalar','real','finite','integer','positive'});
if ~iscell(parameters.amplitudes)||~isequal(size(parameters.amplitudes),[1 2])||...
   ~all(cellfun(@(a)isfloat(a)&&isreal(a)&&all(isfinite(a),'all')&&...
                    isequal(size(a),size(parameters.slice_durs)),parameters.amplitudes))
    error('Spinach:dnp_train:badAmplitudes',...
          'parameters.amplitudes must contain two finite real rows matching parameters.slice_durs.');
end
end


