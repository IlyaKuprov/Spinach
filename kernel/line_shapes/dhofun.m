% Normalised damped harmonic oscillator response function in mag-
% netic resonance notation. This is the standard shape of a phonon
% band: a resonance at the natural frequency of the oscillator, br-
% oadened by damping, and vanishing quadratically at zero frequen-
% cy, which a Lorentzian centred at the same place does not do. In
% the weak damping limit the function tends to a Lorentzian of the
% same width. Evaluation uses frequencies scaled by nat_freq to
% avoid intermediate overflow in single precision. Syntax:
%
%                    y=dhofun(x,nat_freq,fwhm)
%
% Parameters:
%
%            x - argument values, a real array of any dimension;
%                the function is zero at non-positive arguments
%
%     nat_freq - natural frequency of the undamped oscillator, a
%                positive real scalar, in the same units as x; the
%                maximum of the function sits exactly here
%
%         fwhm - damping rate of the oscillator, a positive real
%                scalar, in the same units as x; for this respon-
%                se function it is exactly the full width at half-
%                maximum, at any damping
%
% Outputs:
%
%            y - function values at the points specified in x, an
%                array of the same size and type as x, normalised
%                to unit integral over positive arguments
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=dhofun.m>

function y=dhofun(x,nat_freq,fwhm)

% Check consistency
grumble(x,nat_freq,fwhm);

% The response function only lives at positive frequencies
y=zeros(size(x),'like',x); pos_args=(x>0); arg=x(pos_args);

% Compute the response in reduced frequency units
rel_freq=arg/nat_freq; rel_fwhm=fwhm/nat_freq;
y(pos_args)=(2*rel_fwhm/(pi*nat_freq))*(rel_freq.^2)./...
            ((rel_freq.^2-1).^2+(rel_fwhm*rel_freq).^2);

end

% Consistency enforcement
function grumble(x,nat_freq,fwhm)
if (~isnumeric(x))||(~isreal(x))
    error('x must be an array of real numbers.');
end
if (~isnumeric(nat_freq))||(~isreal(nat_freq))||...
   (numel(nat_freq)~=1)||(~isfinite(nat_freq))||(nat_freq<=0)
    error('nat_freq must be a finite positive real scalar.');
end
if (~isnumeric(fwhm))||(~isreal(fwhm))||...
   (numel(fwhm)~=1)||(~isfinite(fwhm))||(fwhm<=0)
    error('fwhm must be a finite positive real scalar.');
end
end

