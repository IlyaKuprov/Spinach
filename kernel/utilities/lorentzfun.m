% Normalised Lorentzian function in magnetic resonance notation
% with a phase distortion. Syntax:
%
%    [real_part,imag_part]=lorentzfun(offs,ampl,fwhm,x,phi)
%
% Parameters:
%
%     offs - peak offset from zero
%
%     ampl - amplitude multiplier, scalar
%
%     fwhm - full width at half-maximum, scalar
%
%        x - argument, array of any dimension
%
%      phi - phase distortion, radians
%       
% Outputs:
%
%     real_part - an array of values, same size as x
%
%     imag_part - an array of values, same size as x
%              
% u.rasulov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=lorentzfun.m>

function [real_part,imag_part]=lorentzfun(offs,ampl,fwhm,x,phi)

% Check consistency
grumble(offs,ampl,fwhm,x,phi);

% Width parameter
gam=fwhm/2;

% Calculate output
real_part=((ampl/(2*pi*gam))./(1+((x-offs)/gam).^2)).*cos(phi)-...
          ((x-offs)/gam).*((ampl/(2*pi*gam))./(1+((x-offs)/gam).^2)).*sin(phi);
imag_part=((ampl/(2*pi*gam))./(1+((x-offs)/gam).^2)).*sin(phi)+...
          ((x-offs)/gam).*((ampl/(2*pi*gam))./(1+((x-offs)/gam).^2)).*cos(phi);

end

% Consistency enforcement
function grumble(offs,ampl,fwhm,x,phi)
if (~isnumeric(x))||(~isreal(x))
    error('x must be an array of real numbers.');
end
if (~isnumeric(fwhm))||(~isreal(fwhm))||...
   (numel(fwhm)~=1)||(fwhm<=0)
    error('fwhm must be a positive real scalar.');
end
if (~isnumeric(offs))||(~isreal(offs))||(numel(offs)~=1)
    error('offs must be a real scalar.');
end
if (~isnumeric(ampl))||(~isreal(ampl))||(numel(ampl)~=1)
    error('ampl must be a real scalar.');
end
if (~isnumeric(phi))||(~isreal(phi))||(numel(phi)~=1)
    error('ampl must be a real scalar.');
end
end

% Life is so much better if you have a 
% blender - you can blend anything.
%
% Anupama

