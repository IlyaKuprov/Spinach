% Amplifier compression distortion model. Applies a saturating 
% root-sigmoidal distortion to the waveform amplitude:
%
%                     y=x/(1+(x/a)^s)^(1/s)
%
% Treats odd channels of multi-channel waveform as X and even
% ones as Y components; the autodiff Jacobian is returned for
% the vectorisation of the input array. Syntax:
%
%                 [w,J]=amp_root(w,sat_lvls,s)
%
% Parameters:
%
%    w         - waveform in rad/s nutation frequency units,
%                one time slice per column, and rows arran-
%                ged as XYXY... with respect to in-phase and
%                quadrature parts on each control channel
%
%    sat_lvls  - saturation levels beyond which the amplifi-
%                er cannot go, one value per X,Y pair in w,
%                giving the maximum output sqrt(X^2+Y^2)
%
%    s         - a vector of positive integers (one value per
%                X,Y pair in w) regulating the sharpness of
%                the transition from linear to saturating be-
%                haviour, a good starting choice is 4
%
% Outputs:
%
%    w         - distorted waveform in the same units and 
%                layout as the input
%
%    J         - distortion Jacobian matrix with respect to
%                the vectorisation of the input, sparse
%
% u.rasulov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=amp_root.m>

function [w,J]=amp_root(w,sat_lvls,s)

% Check consistency
grumble(w,sat_lvls,s);

% Autodiff wrapper
if nargout<2
    
    % Plain call
    w=distort(w,sat_lvls,s);

else

    % Autodiff call including Jacobian
    [w,J]=dlfeval(@distort,dlarray(w),sat_lvls,s);

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
function [w_dist,J]=distort(w,sat_lvls,s)

% Preallocate output
w_dist=zeros(size(w),'like',w);

% Loop over channels
for n=1:(size(w,1)/2)

    % Get amplitude and phase
    amp=sqrt(w(2*n-1,:).^2+w(2*n,:).^2);
    phi=atan2(w(2*n,:),w(2*n-1,:));

    % Distort the amplitude
    amp=amp./(1+(amp/sat_lvls(n)).^s(n)).^(1/s(n));

    % Get X and Y components back
    w_dist(2*n-1,:)=amp.*cos(phi); 
    w_dist(2*n,:)=amp.*sin(phi);

end

% Compute non-zero blocks of the Jacobian
if nargout>1, J=dljacobian(w_dist,w,1); end

end

% Consistency enforcement
function grumble(w,sat_lvls,s)
if (~isnumeric(w))||(~isreal(w))||(mod(size(w,1),2)~=0)
    error('w must be an array of reals with an even number of rows.');
end
if (~isnumeric(sat_lvls))||(~isreal(sat_lvls))||...
   (numel(sat_lvls)~=size(w,1)/2)||any(sat_lvls<=0,'all')
    error('sat_lvls must be a real array with one element per XY channel pair.');
end
if (~isnumeric(s))||(~isreal(s))||(numel(s)~=size(w,1)/2)||...
   any(mod(s,1)~=0,'all')||any(s<1,'all')
    error('s must be an array of positive integers with one element per XY channel pair.');
end
end

% Mathematicians are like Frenchmen: whatever you say to them
% they translate into their own language and forthwith it is
% something entirely different.
%
% Johann Wolfgang von Goethe

