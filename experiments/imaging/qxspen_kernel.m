% Distortion kernel of the QxSPEN experiment and its derivatives. Syntax:
%
%   [K,dK_dgam,dK_ddel]=qxspen_kernel(FOVy,NSR,Nyacq,alp,bet,gam,del)
%
% Parameters:
%
%     FOVy    - field of view along Y, mm
%
%     NSR     - number of points in the reconstruction
%               by regularisation
%
%     Nyacq   - number of points that is acquired
%               by the instrument
%
%     alp     - an uncertain magic number that IK does
%               not understand, ask Ke Dai
%
%     bet     - an uncertain magic number that IK does
%               not understand, ask Ke Dai
%
%     gam     - slope of the linear phase, rad/mm
%
%     del     - constant phase, radians
%
% Outputs:
%
%     K       - QxSPEN kernel matrix
%
%     dK_dgam - derivative of K with respect to gam
%
%     dK_ddel - derivative of K with respect to del
%
% Note: a dangerous numerical integration stage with a fixed point count
%       is used - must be replaced with an analytical expression!
%
% ke.dai@weizmann.ac.il
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=qxspen_kernel.m>

function [K,dK_dgam,dK_ddel]=qxspen_kernel(FOVy,NSR,Nyacq,alp,bet,gam,del)

% Check consistency
grumble(FOVy,NSR,Nyacq,alp,bet,gam,del);

% Location grid in the [...] dimension
yacq=linspace(-FOVy/2,FOVy/2,Nyacq)';

% Location grid and step in the [...] dimension
ySR=linspace(-FOVy/2,FOVy/2,NSR); deltaySR=ySR(2)-ySR(1);

% Preallocate the answers
K=zeros([Nyacq NSR],'like',1i);

% Loop over [...]
for m=1:NSR

    % Location grid in the integration dimension
    yint=linspace(ySR(m)-deltaySR/2,...
                  ySR(m)+deltaySR/2,100); % TODO: adaptive accuracy control

    % Compute kernel element
    K(:,m)=trapz(sinc(alp*(yint-yacq)).*...         % Sinc part of the integrand
                 exp(1i*bet*(yint-yacq).^2),2).*... % Exp part of the integrand
                 exp(1i*bet*yacq.^2).*...           % Quadratic phase
                 exp(1i*gam*yacq)*...               % Linear phase
                 exp(1i*del);                       % Constant phase

end

% Compute the derivatives
dK_dgam=1i*yacq.*K; dK_ddel=1i*K;

% Normalise the outputs
nfactor=norm(K,2); K=K/nfactor;
dK_dgam=dK_dgam/nfactor;
dK_ddel=dK_ddel/nfactor;

end

% Consistency enforcement
function grumble(FOVy,NSR,Nyacq,alp,bet,gam,del)
if (~isnumeric(FOVy))||(~isreal(FOVy))||(~isscalar(FOVy))||...
   (~isnumeric(NSR))||(~isreal(NSR))||(~isscalar(NSR))||...
   (~isnumeric(Nyacq))||(~isreal(Nyacq))||(~isscalar(Nyacq))||...
   (~isnumeric(alp))||(~isreal(alp))||(~isscalar(alp))||...
   (~isnumeric(bet))||(~isreal(bet))||(~isscalar(bet))||...
   (~isnumeric(gam))||(~isreal(gam))||(~isscalar(gam))||...
   (~isnumeric(del))||(~isreal(del))||(~isscalar(del))
    error('all arguments must be real scalars.');
end
end

% Editors are short-sighted fear-based 
% life forms.
%
% Wednesday Addams

