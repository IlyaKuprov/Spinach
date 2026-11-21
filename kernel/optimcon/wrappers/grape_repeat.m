% Evaluate GRAPE with repeated, tied blocks of Cartesian controls.
% Syntax:
%
%   [traj_data,fidelity,grad,hess]=grape_repeat(waveform,spin_system)
%
% Parameters:
%
%   waveform     - normalised control values, one column per block sample
%   spin_system  - full-duration control problem initialised by optimcon
%
% Outputs:
%
%   traj_data    - full-duration trajectory data
%   fidelity     - physical fidelity and separate full-waveform penalties
%   grad,hess    - derivatives with respect to the short waveform
%
% The block length must divide the full number of waveform columns.
% The full timing grid, power scaling, and penalties are unchanged.
% Use explicit Cartesian samples without a waveform basis or freeze
% mask; independent freezes cannot generally be tied across repeats.
%
% talos@spindynamics.org

function [traj_data,fidelity,grad,hess]=grape_repeat(waveform,spin_system)

% Check consistency and expand the tied controls
grumble(waveform,spin_system);
if nargout<2
    error('request at least trajectory and fidelity outputs.');
end
nreps=spin_system.control.pulse_ntpts/size(waveform,2);
expanded=repmat(waveform,1,nreps);

% Evaluate the full waveform with the requested derivative order
switch nargout
    case 2
        [traj_data,fidelity]=grape_xy(expanded,spin_system);
    case 3
        [traj_data,fidelity,grad]=grape_xy(expanded,spin_system);
    case 4
        [traj_data,fidelity,grad,hess]=grape_xy(expanded,spin_system);
end

% Sum each full-time gradient over the tied copies
if nargout>=3
    grad=reshape(grad,size(waveform,1),size(waveform,2),nreps,[]);
    grad=reshape(sum(grad,3),size(waveform,1),size(waveform,2),[]);
end

% Contract both axes of each full-time Hessian
if nargout==4
    map=repmat(speye(numel(waveform)),nreps,1);
    pulled=zeros(numel(waveform),numel(waveform),size(hess,3));
    for n=1:size(hess,3)
        pulled(:,:,n)=map'*hess(:,:,n)*map;
    end
    hess=pulled;
end

end

% Consistency enforcement
function grumble(waveform,spin_system)
if ~isfield(spin_system,'control')
    error('call optimcon on the full control problem first.');
end
if ~isnumeric(waveform)||~isreal(waveform)||~ismatrix(waveform)||...
   isempty(waveform)||any(~isfinite(waveform),'all')||...
   size(waveform,1)~=spin_system.control.ncontrols
    error('waveform must contain finite real samples, one row per control.');
end
if mod(spin_system.control.pulse_ntpts,size(waveform,2))~=0
    error('the block length must divide the full waveform length.');
end
if ~isempty(spin_system.control.basis)||~isempty(spin_system.control.freeze)
    error('tied repetition requires explicit samples without a freeze mask.');
end
end


