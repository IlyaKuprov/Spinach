% Normalised Gaussian function in magnetic resonance notation and
% its convolution with a triangular function. Syntax:
%
%                    y=gausscon(offs,ampl,fwhm,x)
%
% Parameters:
%
%     offs - peak offset from zero - when this is a scalar,
%            a Gaussian is returned; when this is a vector
%            with three elements, a convolution with a tri-
%            angular function is returned.
%
%     ampl - amplitude multiplier, scalar
%
%     fwhm - full width at half-maximum, scalar
%
%        x - argument, array of any dimension
%
% Outputs:
%
%        y - an array of values, same size as x
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=gausscon.m>

function y=gausscon(offs,ampl,fwhm,x)

% Check consistency
grumble(offs,ampl,fwhm,x);

% Compute standard deviation
sigma=fwhm/(2*sqrt(2*log(2)));

% Sort the offsets
offs=sort(offs(:),1,'ascend');

% Similarity tolerance
sim_tol=sqrt(eps)*norm(offs,2);
if sim_tol==0
    sim_tol=eps;
end

% Single offset or three identical vertices
if isscalar(offs)||(offs(3)-offs(1)<=sim_tol)

    % Just a Gaussian curve
    y=ampl*gaussfun(x-mean(offs),fwhm);

% Vertices 1 and 2 identical
elseif (offs(2)-offs(1)<=sim_tol)

    % Update the offsets
    offs=[mean(offs(1:2)) offs(3)];

    % Compute Gaussian values at the triangle limits
    gauss_a=gaussfun(x-offs(1),fwhm);
    gauss_b=gaussfun(x-offs(2),fwhm);

    % Integrate the Gaussian between the triangle limits
    prob_int=(erf((offs(2)-x)/(sqrt(2)*sigma))-...
              erf((offs(1)-x)/(sqrt(2)*sigma)))/2;

    % Gaussian convolution with a right-angle triangle
    y=2*ampl*((offs(2)-x).*prob_int+sigma^2*(gauss_b-gauss_a))/...
      ((offs(2)-offs(1))^2);

% Vertices 2 and 3 identical
elseif (offs(3)-offs(2)<=sim_tol)

    % Update the offsets
    offs=[offs(1) mean(offs(2:3))];

    % Compute Gaussian values at the triangle limits
    gauss_a=gaussfun(x-offs(1),fwhm);
    gauss_b=gaussfun(x-offs(2),fwhm);

    % Integrate the Gaussian between the triangle limits
    prob_int=(erf((offs(2)-x)/(sqrt(2)*sigma))-...
              erf((offs(1)-x)/(sqrt(2)*sigma)))/2;

    % Gaussian convolution with a right-angle triangle
    y=2*ampl*((x-offs(1)).*prob_int-sigma^2*(gauss_b-gauss_a))/...
      ((offs(2)-offs(1))^2);

else

    % Compute Gaussian values at the triangle vertices
    gauss_a=gaussfun(x-offs(1),fwhm);
    gauss_b=gaussfun(x-offs(2),fwhm);
    gauss_c=gaussfun(x-offs(3),fwhm);

    % Integrate the Gaussian over each triangle segment
    prob_ab=(erf((offs(2)-x)/(sqrt(2)*sigma))-...
             erf((offs(1)-x)/(sqrt(2)*sigma)))/2;
    prob_bc=(erf((offs(3)-x)/(sqrt(2)*sigma))-...
             erf((offs(2)-x)/(sqrt(2)*sigma)))/2;

    % Compute first moments over the triangle segments
    left_int=(x-offs(1)).*prob_ab-sigma^2*(gauss_b-gauss_a);
    right_int=(offs(3)-x).*prob_bc+sigma^2*(gauss_c-gauss_b);

    % Gaussian convolution with a general triangle
    y=2*ampl*(left_int/((offs(2)-offs(1))*(offs(3)-offs(1)))+...
              right_int/((offs(3)-offs(2))*(offs(3)-offs(1))));

end

end

% Consistency enforcement
function grumble(offs,ampl,fwhm,x)
if (~isnumeric(offs))||(~isreal(offs))||...
   (~ismember(numel(offs),[1 3]))||any(~isfinite(offs(:)))
    error('offs must be a finite real scalar or a three-element vector.');
end
if (~isnumeric(ampl))||(~isreal(ampl))||...
   (numel(ampl)~=1)||(~isfinite(ampl))
    error('ampl must be a finite real scalar.');
end
if (~isnumeric(x))||(~isreal(x))||any(~isfinite(x(:)))
    error('x must be an array of finite real numbers.');
end
if (~isnumeric(fwhm))||(~isreal(fwhm))||...
   (numel(fwhm)~=1)||(~isfinite(fwhm))||(fwhm<=0)
    error('fwhm must be a finite positive real number.');
end
end

% As for our majority... one is enough. 
%
% Benjamin Disraeli

