% Applies sampled complex RF gains to interleaved X/Y control rows.
% Syntax:
%
%                    [w,J]=rotor_rf(w,gain)
%
% Parameters:
%
%    w    - real [2*nchannels x nticks] waveform
%    gain - complex [nchannels x nticks] dimensionless gains
%
% Outputs:
%
%    w - modulated waveform with the original dimensions
%    J - sparse real Jacobian acting on w(:) in column-major order
%
% Talos <talos@spindynamics.org>

function [w,J]=rotor_rf(w,gain)

% Check consistency
grumble(w,gain);

% Multiply each complex quadrature pair by its sampled gain
rf=gain.*complex(w(1:2:end,:),w(2:2:end,:));
w(1:2:end,:)=real(rf);
w(2:2:end,:)=imag(rf);

% Assemble one real two-by-two block per complex sample
if nargout>1

    % Use consecutive vector indices for each quadrature pair
    x=(1:2:numel(w)).';
    y=x+1;
    J=sparse([x;y;y;x],[x;y;x;y],...
             double([real(gain(:));real(gain(:));imag(gain(:));-imag(gain(:))]),...
             numel(w),numel(w));

end

end

% Consistency enforcement
function grumble(w,gain)
if ~isfloat(w)||~ismatrix(w)||isempty(w)||~isreal(w)||any(~isfinite(w),'all')
    error('w must be a nonempty finite real floating-point matrix.');
end
if ~isfloat(gain)||~isequal(2*size(gain,1),size(w,1))||...
   ~isequal(size(gain,2),size(w,2))||~ismatrix(gain)||any(~isfinite(gain),'all')
    error('gain must contain one finite complex gain per quadrature pair and tick.');
end
end


