% Amplifier compression distortion model. Applies a saturating 
% hyperbolic tangent distortion to the user-supplied waveform:
%
%                         y=a*tanh(x/a)
%
% Treats odd channels of multi-channel waveform as X and even 
% ones as Y components; the autodiff Jacobian is returned for
% the vectorisation of the input array. Syntax:
%
%                  [w,J]=amp_tanh(w,sat_lvls)
%
% Parameters:
%
%    w         - waveform in rad/s nutation frequency units,
%                one time slice per column, and rows arran-
%                ged as XYXY... with respect to in-phase and
%                quadrature parts on each control channel
%
%    sat_lvls  - saturation level beyond which the amplifier
%                cannot go, rad/s nutation frequency units;
%                one value per XY pair in w, corresponding 
%                to the maximum output sqrt(X^2+Y^2) value
%
% Outputs:
%
%    w         - distorted waveform in the same units and 
%                layout as the input
%
%    J         - distortion Jacobian matrix with respect to
%                the vectorisation of the input, sparse
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=amp_tanh.m>

function [w,J]=amp_tanh(w,sat_lvls)

% Autodiff wrapper
if nargout<2
    
    % Plain call
    w=distort(w,sat_lvls);

else

    % Autodiff call including Jacobian
    [w,J]=dlfeval(@distort,dlarray(w),sat_lvls);

    % Strip the autodiff rigging
    w=extractdata(w); J=extractdata(J);

    % Flatten out the Jacobian
    J=squeeze(num2cell(J,[1 2])); 
    for n=1:numel(J)
        J{n}=sparse(J{n});
    end
    J=blkdiag(J{:});

end

end

% Actual distortion function
function [w_dist,J]=distort(w,sat_lvls)

% Preallocate output
w_dist=zeros(size(w),'like',w);

% Loop over channels
for n=1:(size(w,1)/2)

    % Get amplitude and phase
    amp=sqrt(w(2*n-1,:).^2+w(2*n,:).^2);
    phi=atan2(w(2*n,:),w(2*n-1,:));

    % Distort the amplitude
    amp=sat_lvls(n)*tanh(amp/sat_lvls(n));

    % Get X and Y components back
    w_dist(2*n-1,:)=amp.*cos(phi); 
    w_dist(2*n,:)=amp.*sin(phi);

end

% Compute non-zero blocks of the Jacobian
if nargout>1, J=dljacobian(w_dist,w,1); end

end

% In mathematics you don't understand things. You
% just get used to them.
%
% John von Neumann

