% Chirp pulse waveform with a sine bell power or a quarter-sine 
% amplitude fade-in and fade-out. Generates unidirectional chir-
% ps or saltire chirps which are super-positions of two counter-
% sweeping chirps. Syntax:
%
%    [Cx,Cy,durs,ints,amps,phis,frqs]=...
%                     chirp_pulse(npts,dur,bwidth,smp,type)
%
% Parameters:
%
%        npts    - number of discretization points in
%                  the waveform
%
%    duration    - pulse duration, seconds
%
%      bwidth    - chirp sweep bandwidth around
%                  zero frequency, Hz
%
%        type    - 'wurst', 'smoothed', or 'saltire'; the
%                   default is uniform time grid, to get
%                   adaptive sampling, add '-adaptive'
%
%         smp    - smoothing parameter; for 'wurst', this
%                  is the power in 
%           
%                               1-|sin(x)^smp|
%
%                  as x approaches pi/2 at either the edge
%                  of the pulse. For 'smoothed' and 'salti-
%                  re', this is the fraction of the pulse
%                  duration (in percent) that is affected
%                  by a sine bell fade-in and fade-out: 0
%                  means square amplitude envelope and 50
%                  means sine bell envelope.
%
% Outputs:
%
%          Cx    - real part of the waveform, calibrated to
%                  produce an inversion pulse, rad/s
%
%          Cy    - imag part of the waveform, calibrated to
%                  produce an inversion pulse, rad/s
%
%        durs    - slice durations for piecewise-constant
%                  approximation, seconds
%
%        ints    - interval durations for piecewise-linear
%                  approximation, seconds
%
%        amps    - waveform amplitudes, rad/s
%
%        phis    - waveform phases, rad
%
%        frqs    - waveform frequencies, Hz
%
%   intv_grid    - normalised interval grid, npts-1 elements
%
% Note: Cy is zero for the saltire pulse, this radically changes
%       its phase and amplitude profiles.
%
% ilya.kuprov@weizmann.ac.uk
% jeannicolas.dumez@cnrs.fr
% ludmilla.guduff@cnrs.fr
% m.g.concilio@soton.ac.uk
% mohammadali.foroozandeh@chem.ox.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=chirp_pulse.m>

function [Cx,Cy,durs,ints,amps,phis,frqs]=...
                     chirp_pulse(npts,dur,bwidth,smp,type)

% Check consistency
grumble(npts,dur,bwidth,smp,type);

% Recursive call for saltire
if contains(type,'saltire')

    % Saltire pulses are not monochromatic
    if nargout>6, error('saltire pulse is not monochromatic.'); end

    % Request a smoothed chirp
    type=replace(type,'saltire','smoothed');
    [Cx,Cy,durs,ints]=chirp_pulse(npts,dur,bwidth,smp,type);

    % Zero out Y component and update parameters
    Cy=zeros(size(Cy)); amps=abs(Cx); phis=(pi/2)*sign(Cy); return;

end

% Decide the time grid
if contains(type,'-adaptive')

    % Adaptive normalised time grid
    time_grid=linspace(-1.0,1.0,npts);
    time_grid=sign(time_grid).*abs(time_grid).^0.75;
    time_grid=time_grid/2;

    % For PWL: interval durations 
    ints=dur*diff(time_grid);

    % For PWC: slice durations
    aux_grid=linspace(-1.0,1.0,npts+1);
    aux_grid=sign(aux_grid).*abs(aux_grid).^0.75;
    aux_grid=aux_grid/2; durs=dur*diff(aux_grid);

else

    % Uniform normalised time grid
    time_grid=linspace(-0.5,0.5,npts);

    % For PWL: interval durations 
    ints=dur*diff(time_grid);

    % For PWC: slice durations
    durs=(dur/npts)*ones(1,npts);

end

% Phase and frequency sequence
phis=pi*dur*bwidth*(time_grid.^2);
frqs=bwidth*time_grid;

% Sampling adequacy check
freq_grid=linspace(-bwidth/2,bwidth/2,npts-1);
phi_jumps=2*pi*abs(freq_grid).*ints;
if any(abs(phi_jumps)>pi,'all')&&(nargout<7)
    error('insufficient number of points to sample the pulse.');
end

% Amplitude sequence
switch type
    
    case {'wurst','wurst-adaptive'}
        
        % Pulse amplitude envelope
        amps=1-abs(sin(pi*time_grid).^smp);
       
    case {'smoothed','smoothed-adaptive'}

        % Number of points to fade either side
        np_do_smooth=nnz(time_grid>(50-smp)/100);

        % Number of points to leave intact
        np_no_smooth=npts-2*np_do_smooth;

        % Compute the amplitude envelope
        fade_in=sin(linspace(0,pi/2,np_do_smooth));
        fade_out=fliplr(fade_in);
        amps=[fade_in ones(1,np_no_smooth) fade_out];
        
    otherwise
        
        % Complain and bomb out
        error('unknown pulse type.');
        
end

% Calibrate pulse amplitudes
amps=2*pi*sqrt(bwidth/dur)*amps;

% Transform into Cartesian representation
[Cx,Cy]=polar2cartesian(amps,phis);

end

% Consistency enforcement
function grumble(npts,dur,bwidth,smp,type)
if (~isnumeric(dur))||(~isreal(dur))||...
   (~isfinite(dur))||(numel(dur)~=1)||(dur<=0)
    error('dur must be a positive real number.');
end
if (~isnumeric(bwidth))||(~isreal(bwidth))||...
   (~isfinite(bwidth))||(numel(bwidth)~=1)||(bwidth<=0)
    error('bwidth must be a positive real number.');
end
if (~isnumeric(npts))||(~isreal(npts))||(~isfinite(npts))||...
   (numel(npts)~=1)||(npts<1)||(mod(npts,1)~=0)
    error('npts must be a positive real integer.');
end
if contains(type,'wurst')&&(smp<1)
    error('smp for ''wurst'' pulse must be greater than 1.');
end
if contains(type,'smoothed')&&((smp<0)||(smp>50))
    error('smp for ''smoothed'' must be between 0 and 50.');
end
if contains(type,'saltire')&&((smp<0)||(smp>50))
    error('smp for ''saltire'' must be between 0 and 50.');
end
end

% Wer ein WARUM zum Leben hat, ertragt fast jedes WIE.
%
% Unofficial motto of Southampton
% Chemistry Department

