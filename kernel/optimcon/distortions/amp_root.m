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

% Call the distortion function
if nargout<2

    % Plain call
    w=distort(w,sat_lvls,s);

else

    % Distortion call including Jacobian
    [w,J]=distort(w,sat_lvls,s);

end

end

% Actual distortion function
function [w_dist,J]=distort(w,sat_lvls,s)

% Preallocate output
w_dist=zeros(size(w));
if nargout>1
    rows=zeros(2*numel(w),1);
    cols=zeros(2*numel(w),1);
    vals=zeros(2*numel(w),1);
    jac_idx=0;
end

% Loop over channels
for n=1:(size(w,1)/2)

    % Loop over time points
    for k=1:size(w,2)

        % Get Cartesian components
        x=w(2*n-1,k); y=w(2*n,k);

        % Get radial amplitude
        amp=sqrt(x^2+y^2);

        % Compute radial scale and curvature
        if amp==0
            scale=1; curv=0;
        else
            amp_pwr=(amp/sat_lvls(n))^s(n);
            scale=(1+amp_pwr)^(-1/s(n));
            curv=-(amp^(s(n)-2))/(sat_lvls(n)^s(n))*...
                 (1+amp_pwr)^(-1/s(n)-1);
        end

        % Apply radial compression
        w_dist(2*n-1,k)=scale*x;
        w_dist(2*n,k)=scale*y;

        % Fill the Cartesian Jacobian block
        if nargout>1
            idx_x=2*n-1+size(w,1)*(k-1);
            idx_y=2*n+size(w,1)*(k-1);
            rows(jac_idx+(1:4))=[idx_x; idx_x; idx_y; idx_y];
            cols(jac_idx+(1:4))=[idx_x; idx_y; idx_x; idx_y];
            vals(jac_idx+(1:4))=[scale+curv*x^2; curv*x*y; ...
                                 curv*x*y; scale+curv*y^2];
            jac_idx=jac_idx+4;
        end

    end

end

% Assemble the sparse Jacobian
if nargout>1
    J=sparse(rows,cols,vals,numel(w),numel(w));
end

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
