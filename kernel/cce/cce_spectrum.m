% Fourier transform of a declared stationary correlation. Syntax:
%
%               spectrum=cce_spectrum(times,correlation,parameters)
%
% times is a strictly increasing seconds column; correlation has one row
% per time. parameters.spectrum requires convention='exp(+i*omega*t)',
% sidedness='two-sided'/'one-sided', window='none'/'hann', frequency_grid
% (angular frequencies rad/s), time_domain='full'/'positive', and
% weights='trapezoid'/'rectangle'. Positive-time input is extended by
% C(-t)=conj(C(t)), and therefore asserts stationarity. One-sided output
% requires a real even correlation and doubles only nonzero frequencies.
% Rectangle uses uniform dt at every sample, suitable for DFT Parseval;
% trapezoid uses endpoint half weights. No window-power renormalization.
% Optional spectrum.modulation is a real column on the ORIGINAL time
% mesh: filter=|integral y(t)exp(+i*omega*t)dt|^2, in seconds squared.
%
% Inputs:
%   times - time samples in seconds.
%   correlation - time-by-observable complex correlations.
%   parameters.spectrum - frequency_grid (rad/s), convention=exp(+i*omega*t),
%   sidedness=two-sided|one-sided, window=none|hann,
%   time_domain=full|positive, weights=trapezoid|rectangle; positive-time
%   extension additionally requires stationary=true.
%
% Outputs:
%   spectrum - angular-frequency grid, transformed values and declared
%   convention/quadrature metadata.
%
% See also cce_noise.
%
% ilya.kuprov@weizmann.ac.il

function spectrum=cce_spectrum(times,correlation,parameters)

% Validate the transform convention, dimensional axes and quadrature
grumble(times,correlation,parameters);
s=parameters.spectrum; times=times(:); omega=s.frequency_grid(:);
original_times=times;

% Extend only an explicitly declared stationary positive-time correlation
if strcmp(s.time_domain,'positive')
    times=[-flipud(times(2:end));times];
    correlation=[conj(flipud(correlation(2:end,:)));correlation];
end

% Construct quadrature and an explicitly unnormalized symmetric window
weights=quadrature(times,s.weights); window=ones(size(times));
if strcmp(s.window,'hann')
    window=(1-cos(2*pi*(times-times(1))/(times(end)-times(1))))/2;
end
weighted=weights.*window.*correlation; values=complex(zeros(numel(omega),size(correlation,2)));

% Evaluate directly without hidden FFT shifts or odd-length assumptions
for k=1:numel(omega)
    values(k,:)=exp(1i*omega(k)*times.')*weighted;
end

% Fold only real even spectra into a one-sided spectral density
if strcmp(s.sidedness,'one-sided')
    values(omega>0,:)=2*values(omega>0,:);
end
spectrum=struct('omega',omega,'frequency_hz',omega/(2*pi),'value',values,...
                'convention',s.convention,'sidedness',s.sidedness,...
                'window',window,'time_weights',weights,'times',times,...
                'frequency_weights',quadrature(omega,'trapezoid'),...
                'normalization','C(0)=integral S(omega)domega/(2*pi)',...
                'units','correlation_units*seconds');

% Return the filter transform separately from any noise-model approximation
if isfield(s,'modulation')
    w=quadrature(original_times,s.weights); filter=zeros(numel(omega),1);
    for k=1:numel(omega)
        filter(k)=abs(exp(1i*omega(k)*original_times.')*(w.*s.modulation(:)))^2;
    end
    spectrum.filter=filter; spectrum.filter_units='seconds^2';
end

end

% Generate endpoint-aware weights on uniform or nonuniform axes
function weights=quadrature(x,rule)
if isscalar(x), weights=0; return; end
if strcmp(rule,'rectangle')
    weights=repmat((x(end)-x(1))/(numel(x)-1),size(x));
else
    weights=[x(2)-x(1);x(3:end)-x(1:end-2);x(end)-x(end-1)]/2;
end
end

% Enforce explicit stationary extensions and one-sided quantum restrictions
function grumble(times,correlation,parameters)
if ~isnumeric(times)||~isreal(times)||~isvector(times)||numel(times)<2||...
   any(~isfinite(times))||any(diff(times)<=0)||...
   ~isnumeric(correlation)||size(correlation,1)~=numel(times)||...
   any(~isfinite(correlation),'all')||~ismatrix(correlation)
    error('Finite increasing times and one correlation row per time are required.');
end
if ~isstruct(parameters)||~isfield(parameters,'spectrum')||...
   ~all(isfield(parameters.spectrum,{'convention','sidedness','window',...
                                   'frequency_grid','time_domain','weights'}))
    error('Explicit spectrum convention, sidedness, window, axes and weights are required.');
end
s=parameters.spectrum;
if ~strcmp(s.convention,'exp(+i*omega*t)')||...
   ~ismember(s.sidedness,{'one-sided','two-sided'})||...
   ~ismember(s.window,{'none','hann'})||...
   ~ismember(s.time_domain,{'positive','full'})||...
   ~ismember(s.weights,{'trapezoid','rectangle'})
    error('Unsupported Fourier convention, sidedness, window, domain or quadrature.');
end
w=s.frequency_grid;
if ~isnumeric(w)||~isreal(w)||~isvector(w)||isempty(w)||...
   any(~isfinite(w))||any(diff(w)<=0)
    error('frequency_grid must be strictly increasing finite angular frequencies.');
end
if strcmp(s.time_domain,'positive')
    if ~isfield(s,'stationary')||~islogical(s.stationary)||...
       ~isscalar(s.stationary)||~s.stationary
        error('Positive-time extension requires scalar logical stationary=true.');
    end
    if times(1)~=0, error('Positive-time stationary input must start at zero.'); end
end
if strcmp(s.weights,'rectangle')&&...
   max(abs(diff(times)-mean(diff(times))))>100*eps*max(abs(times))
    error('Rectangle quadrature requires a uniform time grid.');
end
if strcmp(s.sidedness,'one-sided')
    % Compare dimensionless residuals separately for each correlation column
    scale=max(abs(correlation),[],1); scale(scale==0)=1;
    scaled=correlation./scale;
    if any(w<0)||any(max(abs(imag(scaled)),[],1)>100*eps)
        error('One-sided folding requires nonnegative frequencies and a real even correlation.');
    end
    if strcmp(s.time_domain,'full')&&...
       (norm(times(:)+flipud(times(:)))>100*eps*max(abs(times))||...
        any(max(abs(scaled-flipud(scaled)),[],1)>100*eps))
        error('Full-time one-sided input must have symmetric times and an even correlation.');
    end
end
if isfield(s,'modulation')&&(~isnumeric(s.modulation)||~isreal(s.modulation)||...
   ~isvector(s.modulation)||numel(s.modulation)~=numel(times)||any(~isfinite(s.modulation)))
    error('modulation must be a finite real vector on the original time mesh.');
end
end


