% Exports common Spinach EPR/EMR results using sequence parameters.
% Syntax:
%
%        text=jcamp_epr(spin_system,parameters,signal,kind,info)
%
% Parameters:
%
%    spin_system - Spinach system structure
%
%    parameters  - original or returned sequence parameters; see kind
%
%    signal      - native result array or named component structure;
%                  field/ENDOR scans are rows, time/frequency 1D data
%                  are columns, and 2D data have [F2,F1] array order
%
%    kind        - 'field': fieldsweep() row result, with returned
%                  parameters.b_axis row in tesla and mw_freq in Hz
%                  'endor': ENDOR row result, with n_frq row in Hz
%                  'time': 1D/2D pulse signal, dwell 1/sweep in seconds;
%                  dimension count is numel(parameters.npoints) or, for
%                  HYSCORE, numel(parameters.nsteps)
%                  'frequency': processed 1D/2D spectrum; dimension
%                  count is numel(parameters.zerofill); offsets and
%                  sweep widths are in Hz, scalar or [F1,F2]
%
%    info        - metadata structure documented by jcamp_signal;
%                  method and detection are explicit, not inferred
%
% Outputs:
%
%    text - complete JCAMP character row, also written when requested
%
% Frequencies follow ft_axis and array lengths, not plot-unit conversions.
% Acquisition times are relative to the first simulated sample. Scanned
% fields and RF frequencies are used exactly as returned by Spinach;
% their row storage is explicitly translated into JCAMP column traces.
% ENDOR n_frq is the applied RF coordinate; no absolute-value folding is
% performed. Use jcamp_signal for explicit DEER/ESEEM delays, arbitrary
% echo trajectories, non-uniform sampling, or other array layouts.
%
% talos@spindynamics.org

function text=jcamp_epr(spin_system,parameters,signal,kind,info)

% Check consistency
grumble(parameters,signal,kind);

% Use the coordinate arrays returned by swept experiments
if any(strcmp(kind,{'field','endor'}))
    if strcmp(kind,'field')
        axes={parameters.b_axis.'}; units={'TESLA'}; names={'Magnetic field'};
        info.metadata=[{'.MICROWAVE FREQUENCY',parameters.mw_freq}; info.metadata];
    else
        axes={parameters.n_frq.'}; units={'HERTZ'}; names={'RF frequency'};
    end
    if isstruct(signal)
        channels=fieldnames(signal);
        for n=1:numel(channels), signal.(channels{n})=signal.(channels{n}).'; end
    else
        signal=signal.';
    end
else

    % Determine the declared acquisition or processing dimensionality
    if strcmp(kind,'frequency')
        nd=numel(parameters.zerofill);
    elseif isfield(parameters,'npoints')
        nd=numel(parameters.npoints);
    else
        nd=numel(parameters.nsteps);
    end
    if isstruct(signal)
        components=struct2cell(signal); sizes=size(components{1});
    else
        sizes=size(signal);
    end
    sweep=parameters.sweep;
    if isscalar(sweep), sweep=repmat(sweep,1,nd); end
    order=1:nd;
    if nd==2, order=[2 1]; end
    axes=cell(1,nd); units=cell(1,nd); names=cell(1,nd);
    for n=1:nd
        dimension=order(n);
        if strcmp(kind,'time')
            axes{n}=(0:sizes(n)-1)'/sweep(dimension); units{n}='SECOND';
        else
            offset=parameters.offset;
            if isscalar(offset), offset=repmat(offset,1,nd); end
            axes{n}=ft_axis(offset(dimension),sweep(dimension),sizes(n)).';
            units{n}='HERTZ';
        end
        names{n}=['F' num2str(dimension) ' ' kind];
    end
end

% Build the internal JCAMP structure and call the final writer
text=jcamp_signal(spin_system,axes,signal,units,names,info);

end

% Consistency enforcement
function grumble(parameters,signal,kind)
if ~isstruct(parameters)||~isscalar(parameters)||...
   ~ischar(kind)||~isrow(kind)||~any(strcmp(kind,{'field','endor','time','frequency'}))
    error('parameters must be scalar and kind field, endor, time, or frequency.');
end
if isstruct(signal)
    if ~isscalar(signal)||isempty(fieldnames(signal))
        error('signal must be a non-empty scalar component structure.');
    end
    signals=struct2cell(signal);
else
    signals={signal};
end
if any(strcmp(kind,{'field','endor'}))
    if strcmp(kind,'field'), axis=parameters.b_axis; else, axis=parameters.n_frq; end
    if ~isfloat(axis)||~isreal(axis)||~isrow(axis)||isempty(axis)||any(~isfinite(axis))
        error('the scan coordinate must be a non-empty finite real floating-point row.');
    end
    for n=1:numel(signals)
        if ~isfloat(signals{n})||~isrow(signals{n})||numel(signals{n})~=numel(axis)
            error('scan signals must be floating-point rows matching the scan coordinate.');
        end
    end
else
    if strcmp(kind,'frequency')
        nd=numel(parameters.zerofill);
    elseif isfield(parameters,'npoints')
        nd=numel(parameters.npoints);
    else
        nd=numel(parameters.nsteps);
    end
    if ~ismember(nd,[1 2])
        error('time/frequency EPR arrays must have one or two physical dimensions.');
    end
    sweep=parameters.sweep;
    if ~isnumeric(sweep)||~isreal(sweep)||~isvector(sweep)||...
       ~ismember(numel(sweep),[1 nd])||any(~isfinite(sweep))||any(sweep<=0)
        error('sweep must have one positive finite real width per dimension, or be scalar.');
    end
    if strcmp(kind,'frequency')
        offset=parameters.offset;
        if ~isnumeric(offset)||~isreal(offset)||~isvector(offset)||...
           ~ismember(numel(offset),[1 nd])||any(~isfinite(offset))
            error('offset must have one finite real value per dimension, or be scalar.');
        end
    end
    for n=1:numel(signals)
        if ~isfloat(signals{n})||isempty(signals{n})
            error('signal components must be non-empty floating-point arrays.');
        end
        if strcmp(kind,'frequency')&&any(size(signals{n},1:nd)<=2)
            error('frequency dimensions require at least three points for ft_axis.');
        end
    end
end
end


