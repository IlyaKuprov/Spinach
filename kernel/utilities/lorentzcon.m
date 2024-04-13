% Normalised Lorentzian function in magnetic resonance notation and 
% its convolution with a triangular function. Syntax:
%
%                   y=lorentzcon(offs,ampl,fwhm,x)
%
% Parameters:
%
%     offs - peak offset from zero - when this is a scalar,
%            a Lorentzian is returned; when this is a vector
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
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=lorentzcon.m>

function y=lorentzcon(offs,ampl,fwhm,x)

% Check consistency
grumble(x,fwhm);

% Width parameter
gam=fwhm/2;

% Sort the offsets
offs=sort(offs(:),1,'ascend');

% Similarity tolerance
sim_tol=sqrt(eps)*norm(offs,2);

% Single offset or three identical vertices
if isscalar(offs)||(offs(3)-offs(1)<sim_tol)
    
    % Just a Lorentzian curve
    y=(ampl/(2*pi*gam))./(1+((x-mean(offs))/gam).^2);
    
% Vertices 1 and 2 identical
elseif (offs(2)-offs(1)<sim_tol)
    
    % Update the offsets
    offs=[mean(offs(1:2)) offs(3)];
    
    % Lorentzian convolution with a right-angle triangle
    y=(ampl/pi)*((x-offs(2))/((offs(1)-offs(2))*(offs(1)-offs(2)))).*atan2(offs(1)-x,gam)+...
      (ampl/pi)*((offs(2)-x)/((offs(1)-offs(2))*(offs(1)-offs(2)))).*atan2(offs(2)-x,gam)+...
      ((ampl*gam/(2*pi))/((offs(1)-offs(2))*(offs(1)-offs(2))))*log(((offs(1)-x).^2+gam^2)./((offs(2)-x).^2+gam^2));
    
% Vertices 2 and 3 identical
elseif (offs(3)-offs(2)<sim_tol)
    
    % Update the offsets
    offs=[offs(1) mean(offs(2:3))];
    
    % Lorentzian convolution with a right-angle triangle
    y=(ampl/pi)*((offs(1)-x)/((offs(1)-offs(2))*(offs(1)-offs(2)))).*atan2(offs(1)-x,gam)+...
      (ampl/pi)*((x-offs(1))/((offs(1)-offs(2))*(offs(1)-offs(2)))).*atan2(offs(2)-x,gam)+...
      ((ampl*gam/(2*pi))/((offs(1)-offs(2))*(offs(1)-offs(2))))*log(((offs(2)-x).^2+gam^2)./((offs(1)-x).^2+gam^2));
 
else

    % Lorentzian convolution with a general triangle
    y=(ampl/pi)*((offs(1)-x)/((offs(1)-offs(2))*(offs(1)-offs(3)))).*atan2(offs(1)-x,gam)+...
      (ampl/pi)*((x-offs(2))/((offs(1)-offs(2))*(offs(2)-offs(3)))).*atan2(offs(2)-x,gam)+...
      (ampl/pi)*((x-offs(3))/((offs(1)-offs(3))*(offs(3)-offs(2)))).*atan2(offs(3)-x,gam)+...
      ((ampl*gam/(2*pi))/((offs(1)-offs(2))*(offs(1)-offs(3))))*log(((offs(2)-x).^2+gam^2)./((offs(1)-x).^2+gam^2))+...
      ((ampl*gam/(2*pi))/((offs(1)-offs(3))*(offs(3)-offs(2))))*log(((offs(3)-x).^2+gam^2)./((offs(2)-x).^2+gam^2));

end

end

% Consistency enforcement
function grumble(x,fwhm)
if (~isnumeric(x))||(~isreal(x))
    error('x must be an array of real numbers.');
end
if (~isnumeric(fwhm))||(~isreal(fwhm))||...
   (numel(fwhm)~=1)||(fwhm<=0)
    error('fwhm must be a positive real number.');
end
end

% [...] the sexual revolution was the mockery of the century,
% because now women were giving to men for free what they used
% to have to marry us for.
%
% A feminist web site

