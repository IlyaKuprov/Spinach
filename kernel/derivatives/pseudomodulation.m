% Pseudomodulation of uniformly sampled spectra using the Hyde
% et al. Fourier-domain algorithm. Syntax:
%
%       output=pseudomodulation(field,spectrum,mod_amp,mod_order)
%
% Parameters:
%
%       field     - N-by-1 real, ordered, uniformly spaced field
%                   axis
%
%       spectrum  - N-by-M spectrum matrix; rows are field samples,
%                   and columns are independent spectra
%
%       mod_amp   - non-negative modulation amplitude in field units
%
%       mod_order - modulation harmonic order: 0, 1, or 2
%
% Outputs:
%
%       output    - N-by-M pseudomodulated spectrum matrix
%
% The implementation follows Eqs. 5-7 of Hyde et al., J. Magn.
% Reson. 96, 1-13 (1992). After phase-sensitive detection, the
% time-dependent prefactors are set to unity, leaving amplitude
% factors 2i for the first harmonic, and 2 for the second harmonic.
%
% alexey.bogdanov@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=pseudomodulation.m>

function output=pseudomodulation(field,spectrum,mod_amp,mod_order)

% Check consistency
grumble(field,spectrum,mod_amp,mod_order);

% Get field grid parameters
npts=numel(field);
field_step=field(2)-field(1);

% Build the angular frequency axis in Matlab FFT ordering
ang_freq=sign(field_step)*real(fftdiff(1,npts,abs(field_step))/1i);
ang_freq=ang_freq(:);

% Convert the modulation amplitude into Bessel arguments
bessel_arg=mod_amp*ang_freq/2;

% Fourier transform along the field axis
spec_ft=fft(spectrum,npts,1);

% Apply the requested pseudomodulation harmonic
switch mod_order

    case 0

        % Zeroth harmonic
        output=ifft(spec_ft.*besselj(0,bessel_arg),npts,1);

    case 1

        % First harmonic
        output=2i*ifft(spec_ft.*besselj(1,bessel_arg),npts,1);

    case 2

        % Second harmonic
        output=2*ifft(spec_ft.*besselj(2,bessel_arg),npts,1);

end

% Remove round-off imaginary parts from real input spectra
if isreal(spectrum)
    output=real(output);
end

end

% Consistency enforcement
function grumble(field,spectrum,mod_amp,mod_order)
if (~isfloat(field))||isempty(field)||issparse(field)||(~iscolumn(field))
    error('field must be a non-empty dense floating-point column vector.');
end
if (~isreal(field))||any(~isfinite(field(:)))
    error('field must contain finite real numbers.');
end
npts=numel(field);
if npts<3
    error('field must contain at least three points.');
end
field_ref=linspace(field(1),field(end),npts).';
grid_tol=npts*eps(max([1; abs(field)]));
if any(abs(field-field_ref)>grid_tol)
    error('field must be ordered and uniformly spaced.');
end
if field(1)==field(end)
    error('field must have a non-zero sweep width.');
end
if (~isfloat(spectrum))||isempty(spectrum)||issparse(spectrum)||(~ismatrix(spectrum))
    error('spectrum must be a non-empty dense floating-point matrix.');
end
if any(~isfinite(spectrum(:)))
    error('spectrum must not contain NaN or Inf values.');
end
if size(spectrum,1)~=npts
    error('spectrum must have the same number of rows as field has points.');
end
if (~isnumeric(mod_amp))||(~isreal(mod_amp))||(numel(mod_amp)~=1)||...
   (~isfinite(mod_amp))||(mod_amp<0)
    error('mod_amp must be a non-negative finite real scalar.');
end
if (~isnumeric(mod_order))||(~isreal(mod_order))||(numel(mod_order)~=1)||...
   (~ismember(mod_order,[0 1 2]))
    error('mod_order must be 0, 1, or 2.');
end
end
