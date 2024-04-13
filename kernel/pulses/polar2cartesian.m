% Converts [RF_amplitude, RF_phase] representation of a pulse waveform and
% the derivatives of any function with respect to those amplitudes and pha-
% ses into the [RF_x, RF_y] representation and the derivatives of the func-
% tion with respect to those X and Y RF values. Syntax:
%
%  [x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy]=polar2cartesian(r,p,Dr,Dp,Drr,Drp,Dpr,Dpp)
%
% Parameters:
%
%    r    - vector of waveform amplitudes
%
%    p    - vector of waveform phases
%
%    Dr   - optional vector of derivatives of some scalar function
%           with respect to the waveform amplitudes.
%
%    Dp   - optional vector of derivatives of some scalar function
%           with respect to the waveform phases.
%
%    Drr  - matrix of second derivatives of the function with respect
%           to the waveform amplitudes.
%
%    Drp  - matrix of second derivatives of the function with respect
%           to the waveform amplitudes and phases.
%
%    Dpr  - matrix of second derivatives of the function with respect
%           to the waveform phases and amplitudes.
%
%    Dpp  - matrix of second derivatives of the function with respect
%           to the waveform phases.
%
% Outputs:
%
%    x    - vector of waveform amplitudes along X
%
%    y    - vector of waveform amplitudes along Y
%
%    Dx   - vector of derivatives of the function with respect to
%           the waveform amplitudes along X
%
%    Dy   - vector of derivatives of the function with respect to
%           the waveform amplitudes along Y
%
%    Dxx  - optional matrix of second derivatives of a scalar function
%           with respect to the waveform amplitudes along X
%
%    Dxy  - optional matrix of second derivatives of a scalar function
%           with respect to the waveform amplitudes along X and Y
%
%    Dyx  - optional matrix of second derivatives of a scalar function
%           with respect to the waveform amplitudes along Y and X
%
%    Dyy  - optional matrix of second derivatives of a scalar function
%           with respect to the waveform amplitudes along Y
%
% i.kuprov@soton.ac.uk
% david.goodwin@inano.au.dk
%
% <https://spindynamics.org/wiki/index.php?title=polar2cartesian.m>

function [x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy]=polar2cartesian(r,p,Dr,Dp,Drr,Drp,Dpr,Dpp)

% Wrap phase
p=wrapToPi(p);

% Check consistency
if nargin==2
    grumble(r,p);
elseif nargin==4
    grumble(r,p,Dr,Dp);
elseif nargin==8
    grumble(r,p,Dr,Dp,Drr,Drp,Dpr,Dpp);
else
    error('incorrect number of arguments.');
end

% Transform coordinates
x=r.*cos(p); y=r.*sin(p);

% Transform derivatives
if (nargin>2)&&(nargout>2)
    
    Dx=+Dr.*cos(p)-Dp.*sin(p)./r;
    Dy=+Dr.*sin(p)+Dp.*cos(p)./r;
    
end

% Transform second derivatives
if (nargin>4)&&(nargout>4)
    
    Dxx=+diag(((sin(p).^2)./r).*Dr)+diag((sin(2*p)./(r.^2)).*Dp)...
        +(cos(p')*cos(p)).*Drr+((sin(p)./r)'*(sin(p)./r)).*Dpp...
        -(cos(p')*(sin(p)./r)).*Dpr-((sin(p)./r)'*cos(p)).*Drp;
    
    Dxy=-diag((sin(2*p)./(2*r)).*Dr)+diag(((1-2*(cos(p).^2))./(r.^2)).*Dp)...
        +(sin(p')*cos(p)).*Drr-((cos(p)./r)'*(sin(p)./r)).*Dpp...
        -(sin(p')*(sin(p)./r)).*Dpr+((cos(p)./r)'*cos(p)).*Drp;
    
    Dyx=-diag((sin(2*p)./(2*r)).*Dr)-diag(((1-2*(sin(p).^2))./(r.^2)).*Dp)...
        +(cos(p')*sin(p)).*Drr-((sin(p)./r)'*(cos(p)./r)).*Dpp...
        +(cos(p')*(cos(p)./r)).*Dpr-((sin(p)./r)'*sin(p)).*Drp;
    
    Dyy=+diag(((cos(p).^2)./r).*Dr)-diag((sin(2*p)./(r.^2)).*Dp)...
        +(sin(p')*sin(p)).*Drr+((cos(p)./r)'*(cos(p)./r)).*Dpp...
        +(sin(p')*(cos(p)./r)).*Dpr+((cos(p)./r)'*sin(p)).*Drp;
    
end

end

% Consistency enforcement
function grumble(r,p,df_dr,df_dp,d2f_dr2,d2f_drdp,d2f_dpdr,d2f_dp2)
if nargin==2
    if (~isnumeric(r))||(~isreal(r))||(~all(r>=0))
        error('amplitude parameter must be a vector of non-negative real numbers.');
    end
    if (~isnumeric(p))||(~isreal(p))
        error('phase parameter must be a vector of real numbers.');
    end
    if ~all(size(r)==size(p))
        error('amplitude and phase vectors must have the same dimension.');
    end
elseif nargin==4
    if (~isnumeric(r))||(~isreal(r))||(~all(r>=0))
        error('amplitude parameter must be a vector of non-negative real numbers.');
    end
    if (~isnumeric(p))||(~isreal(p))
        error('phase parameter must be a vector of real numbers.');
    end
    if (~isnumeric(df_dr))||(~isreal(df_dr))
        error('df_dA parameter must be a vector of real numbers.');
    end
    if (~isnumeric(df_dp))||(~isreal(df_dp))
        error('df_dphi parameter must be a vector of real numbers.');
    end
    if (~all(size(r)==size(p)))||(~all(size(p)==size(df_dr)))||...
       (~all(size(df_dr)==size(df_dp)))
        error('all input vectors must have the same dimension.');
    end
elseif nargin==7
    if (~isnumeric(r))||(~isreal(r))||(~all(r>=0))
        error('amplitude parameter must be a vector of non-negative real numbers.');
    end
    if (~isnumeric(p))||(~isreal(p))
        error('phase parameter must be a vector of real numbers.');
    end
    if (~isnumeric(df_dr))||(~isreal(df_dr))
        error('df_dA parameter must be a vector of real numbers.');
    end
    if (~isnumeric(df_dp))||(~isreal(df_dp))
        error('df_dphi parameter must be a vector of real numbers.');
    end
    if (~all(size(r)==size(p)))||(~all(size(p)==size(df_dr)))||...
       (~all(size(df_dr)==size(df_dp)))
        error('all input vectors must have the same dimension.');
    end
    if (~isnumeric(d2f_dr2))||(~isreal(d2f_dr2))
        error('d2f_dA2 parameter must be a matrix of real numbers.');
    end
    if (~isnumeric(d2f_dp2))||(~isreal(d2f_dp2))
        error('d2f_dphi2 parameter must be a matrix of real numbers.');
    end
    if (~isnumeric(d2f_drdp))||(~isreal(d2f_drdp))
        error('ddf_dAdphi parameter must be a matrix of real numbers.');
    end
    if (~isnumeric(d2f_dpdr))||(~isreal(d2f_dpdr))
        error('ddf_dAdphi parameter must be a matrix of real numbers.');
    end
    if (size(d2f_dr2,2)~=length(df_dr))||(size(d2f_dr2,1)~=size(d2f_dr2,2))||...
        all(size(d2f_dr2)~=size(d2f_dp2))||all(size(d2f_dp2)~=size(d2f_drdp))||...
        all(size(d2f_drdp)~=size(d2f_dpdr))
        error('all input matrices must have the same, square dimensions.');
    end
end
end

% A public opinion poll is no substitute for thought.
%
% Warren Buffett

