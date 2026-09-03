% Converts the [RF_x, RF_y] representation of a pulse waveform and the
% derivatives of any function with respect to those RF values into the
% [RF_amplitude, RF_phase] representation and the derivatives of the
% function with respect to amplitudes and phases. Syntax:
%
%    [r,p,Dr,Dp,Drr,Drp,Dpr,Dpp]=...
%                       cartesian2polar(x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy)
%
% Parameters:
%
%    x    - vector of waveform amplitudes along X
%
%    y    - vector of waveform amplitudes along Y
%
%    Dx   - optional vector of derivatives of a scalar function
%           with respect to the waveform amplitudes along X
%
%    Dy   - optional vector of derivatives of a scalar function
%           with respect to the waveform amplitudes along Y
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
% Outputs:
%
%    r    - vector of waveform amplitudes
%
%    p    - vector of waveform phases
%
%    Dr   - vector of derivatives of the function with respect
%           to the waveform amplitudes.
%
%    Dp   - vector of derivatives of the function with respect
%           to the waveform phases.
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
% ilya.kuprov@weizmann.ac.il
% david.goodwin@inano.au.dk
%
% <https://spindynamics.org/wiki/index.php?title=cartesian2polar.m>

function [r,p,Dr,Dp,Drr,Drp,Dpr,Dpp]=cartesian2polar(x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy)

% Check consistency
if nargin==2
    grumble(x,y);
elseif nargin==4
    grumble(x,y,Dx,Dy);
elseif nargin==8
    grumble(x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy);
else
    error('incorrect number of arguments.');
end

% Transform coordinates
r=sqrt(x.^2+y.^2); p=atan2(y,x);

% Transform derivatives
if (nargin>2)&&(nargout>2)
    
    % Radius
    Dr=+cos(p).*Dx +sin(p).*Dy;
    
    % Phase
    Dp=-y.*Dx +x.*Dy;
    
end

% Transform second derivatives
if (nargin>4)&&(nargout>4)

    % Row copies for the outer products below
    pr=reshape(p,1,[]); xr=reshape(x,1,[]); yr=reshape(y,1,[]);
    Dxr=reshape(Dx,1,[]); Dyr=reshape(Dy,1,[]);

    % Radius, radius
    Drr=+(cos(pr')*cos(pr)).*Dxx +(sin(pr')*sin(pr)).*Dyy...
        +(sin(pr')*cos(pr)).*Dyx +(cos(pr')*sin(pr)).*Dxy;
    
    % Radius, phase
    Drp=-diag(sin(pr).*Dxr) +diag(cos(pr).*Dyr)...
        -(cos(pr')*yr).*Dxx +(sin(pr')*xr).*Dyy...
        -(sin(pr')*yr).*Dyx +(cos(pr')*xr).*Dxy;
    
    % Phase, radius
    Dpr=-diag(sin(pr).*Dxr) +diag(cos(pr).*Dyr)...
        -(yr'*cos(pr)).*Dxx +(xr'*sin(pr)).*Dyy...
        +(xr'*cos(pr)).*Dyx -(yr'*sin(pr)).*Dxy;
    
    % Phase, phase
    Dpp=-diag(xr.*Dxr) -diag(yr.*Dyr)...
        +(yr'*yr).*Dxx +(xr'*xr).*Dyy...
        -(xr'*yr).*Dyx -(yr'*xr).*Dxy;
    
end

end

% Consistency enforcement
function grumble(x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy)
if nargin==2
    if (~isnumeric(x))||(~isreal(x))
        error('x parameter must be a vector of real numbers.');
    end
    if (~isnumeric(y))||(~isreal(y))
        error('y parameter must be a vector of real numbers.');
    end
    if ~all(size(x)==size(y))
        error('x and y vectors must have the same dimension.');
    end
elseif nargin==4
    if (~isnumeric(x))||(~isreal(x))
        error('x parameter must be a vector of real numbers.');
    end
    if (~isnumeric(y))||(~isreal(y))
        error('y parameter must be a vector of real numbers.');
    end
    if (~isnumeric(Dx))||(~isreal(Dx))
        error('df_dx parameter must be a vector of real numbers.');
    end
    if (~isnumeric(Dy))||(~isreal(Dy))
        error('df_dy parameter must be a vector of real numbers.');
    end
    if (~all(size(x)==size(y)))||(~all(size(y)==size(Dx)))||...
       (~all(size(Dx)==size(Dy)))
        error('all input vectors must have the same dimension.');
    end
elseif nargin==8
    if (~isnumeric(x))||(~isreal(x))
        error('x parameter must be a vector of real numbers.');
    end
    if (~isnumeric(y))||(~isreal(y))
        error('y parameter must be a vector of real numbers.');
    end
    if (~isnumeric(Dx))||(~isreal(Dx))
        error('df_dx parameter must be a vector of real numbers.');
    end
    if (~isnumeric(Dy))||(~isreal(Dy))
        error('df_dy parameter must be a vector of real numbers.');
    end
    if (~isnumeric(Dxx))||(~isreal(Dxx))
        error('d2f_dx2 parameter must be a matrix of real numbers.');
    end
    if (~isnumeric(Dyy))||(~isreal(Dyy))
        error('d2f_dy2 parameter must be a matrix of real numbers.');
    end
    if (~isnumeric(Dxy))||(~isreal(Dxy))
        error('ddf_dxdy parameter must be a matrix of real numbers.');
    end
    if (~isnumeric(Dyx))||(~isreal(Dyx))
        error('ddf_dxdy parameter must be a matrix of real numbers.');
    end
    if (~all(size(x)==size(y)))||(~all(size(y)==size(Dx)))||...
       (~all(size(Dx)==size(Dy)))
        error('all input vectors must have the same dimension.');
    end
    if (size(Dxx,2)~=length(Dx))||(size(Dxx,1)~=size(Dxx,2))||...
        (~isequal(size(Dxx),size(Dyy)))||...
        (~isequal(size(Dyy),size(Dxy)))||...
        (~isequal(size(Dxy),size(Dyx)))
        error('all input matrices must have the same, square dimensions.');
    end
end
end

% Beware of geeks bearing formulas.
%
% Warren Buffett

