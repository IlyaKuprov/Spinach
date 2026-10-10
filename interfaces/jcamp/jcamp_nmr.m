% Exports Spinach 1D, 2D, and 3D NMR results to JCAMP. Syntax:
%
%       text=jcamp_nmr(spin_system,parameters,signal,domains,info)
%
% Parameters:
%
%    spin_system - Spinach structure from create(), containing the field
%
%    parameters  - sequence structure with sweep (Hz), offset (Hz), and
%                  spins. Scalar specifications apply to all dimensions,
%                  as in homonuclear 2D examples; otherwise one entry
%                  per physical dimension F1, F2, F3 is required
%
%    signal      - numeric FID/spectrum, or a scalar component structure
%                  (pos/neg, cos/sin, pos_pos/pos_neg/neg_pos/neg_neg).
%                  1D: column; 2D: [F2,F1]; 3D: [F1,F2,F3], matching
%                  plot_1d, plot_2d, and plot_3d and their example sets
%
%    domains     - cell row in physical F1,F2,F3 order, each entry 'time'
%                  or 'frequency'. Its length specifies dimensionality,
%                  including singleton dimensions and mixed-domain arrays
%
%    info        - scalar structure with title, origin, owner, filename,
%                  sequence (ASCII pulse-sequence description), and
%                  metadata (N-by-2 cell table; cell(0,2) if empty).
%                  filename='' returns text without writing a file.
%                  Any time dimension also requires delay: [RD ID] in
%                  microseconds, and acquisition: 'SIMULTANEOUS',
%                  'SEQUENTIAL', or 'SINGLE'. These are not guessed
%
% Outputs:
%
%    text - complete JCAMP character row, also written when requested
%
% Time coordinates start at acquisition time zero with dwell 1/sweep;
% supplied delays describe acquisition timing, not an axis translation.
% Frequency coordinates use ft_axis, matching Spinach spectral plots;
% signal sizes determine lengths, including zero filling. Hz are retained
% irrespective of parameters.axis_units. Per-dimension nuclei, observation
% frequencies (MHz), offsets (Hz), and domains are recorded in private
% metadata. No FFT, processing, or quadrature recombination is performed.
% Mixed-domain blocks use the type of their tabulated (first array) axis.
% Use jcamp_grid and jcamp_export for non-standard NMR sampling/layouts.
%
% talos@spindynamics.org

function text=jcamp_nmr(spin_system,parameters,signal,domains,info)

% Check consistency
grumble(spin_system,parameters,signal,domains,info);

% Map physical dimensions to the established Spinach array layout
nd=numel(domains); order=1:nd;
if nd==2, order=[2 1]; end
if isstruct(signal)
    components=struct2cell(signal); sizes=size(components{1});
else
    sizes=size(signal);
end
sizes(end+1:nd)=1;
sweep=parameters.sweep; offset=parameters.offset; nuclei=parameters.spins;
if isscalar(sweep), sweep=repmat(sweep,1,nd); end
if isscalar(offset), offset=repmat(offset,1,nd); end
if isscalar(nuclei), nuclei=repmat(nuclei,1,nd); end

% Construct axes in array order without altering the simulation result
axes=cell(1,nd); units=cell(1,nd); names=cell(1,nd);
for n=1:nd
    dimension=order(n);
    if strcmp(domains{dimension},'time')
        axes{n}=(0:sizes(n)-1)'/sweep(dimension); units{n}='SECONDS';
    else
        axes{n}=ft_axis(offset(dimension),sweep(dimension),sizes(n)).';
        units{n}='HZ';
    end
    names{n}=['F' num2str(dimension) ' ' nuclei{dimension} ' ' domains{dimension}];
end

% Derive nuclear observation frequencies from the actual static field
observe=zeros(1,nd);
for n=1:nd
    observe(n)=abs(spin(nuclei{n})*spin_system.inter.magnet/(2*pi))*1e-6;
end
metadata={'.OBSERVE FREQUENCY',observe(nd); '.OBSERVE NUCLEUS',['^' nuclei{nd}];...
          '.PULSE SEQUENCE',info.sequence; '$AXIS NUCLEI',strjoin(nuclei,', ');...
          '$AXIS OBSERVE MHZ',observe; '$AXIS OFFSETS HZ',offset;...
          '$AXIS DOMAINS',strjoin(domains,', ')};
if strcmp(domains{order(1)},'time')
    data_type='NMR FID';
else
    data_type='NMR SPECTRUM';
end
if any(strcmp(domains,'time'))
    metadata=[metadata; {'.DELAY',info.delay; '.ACQUISITION MODE',info.acquisition}];
end

% Build the writer input and delegate the final serialisation
block=struct('title',info.title,'type',data_type,'metadata',{[metadata; info.metadata]});
data=struct('title',info.title,'origin',info.origin,'owner',info.owner,...
            'blocks',{jcamp_grid(axes,signal,units,names,block)});
if ~isempty(info.filename), data.filename=info.filename; end
text=jcamp_export(data);

end

% Consistency enforcement
function grumble(spin_system,parameters,signal,domains,info)
if ~isstruct(spin_system)||~isscalar(spin_system)||...
   ~isstruct(parameters)||~isscalar(parameters)
    error('spin_system and parameters must be scalar Spinach structures.');
end
if ~iscell(domains)||~isrow(domains)||isempty(domains)||numel(domains)>3||...
   ~all(cellfun(@(x) ischar(x)&&isrow(x)&&any(strcmp(x,{'time','frequency'})),domains))
    error('domains must contain one to three time/frequency character rows.');
end
if ~isstruct(info)||~isscalar(info)||...
   ~all(isfield(info,{'title','origin','owner','filename','sequence','metadata'}))
    error('info must contain title, origin, owner, filename, sequence, and metadata.');
end
if ~ischar(info.filename)||(~isempty(info.filename)&&~isrow(info.filename))||...
   ~iscell(info.metadata)||size(info.metadata,2)~=2
    error('info.filename must be a character row and metadata an N-by-2 cell table.');
end
if ~ischar(info.sequence)||~isrow(info.sequence)||isempty(strtrim(info.sequence))
    error('info.sequence must be a non-empty pulse-sequence character row.');
end
if any(strcmp(domains,'time'))&&~all(isfield(info,{'delay','acquisition'}))
    error('time-domain export requires info.delay and info.acquisition.');
end
for field={'sweep','offset'}
    values=parameters.(field{1});
    if ~isnumeric(values)||~isreal(values)||~isvector(values)||...
       any(~isfinite(values))||~ismember(numel(values),[1 numel(domains)])
        error('sweep and offset must have one finite real entry per physical dimension, or be scalar.');
    end
end
if any(parameters.sweep<=0)
    error('sweep widths must be positive.');
end
if ~iscell(parameters.spins)||~isvector(parameters.spins)||...
   ~ismember(numel(parameters.spins),[1 numel(domains)])||...
   ~all(cellfun(@(x) ischar(x)&&isrow(x)&&~isempty(regexp(x,'^\d+[A-Z][a-z]?$', 'once')),parameters.spins))
    error('spins must contain nuclear isotope names, one per dimension or one shared nucleus.');
end
if isstruct(signal)
    if ~isscalar(signal)||isempty(fieldnames(signal))
        error('signal must be a non-empty scalar component structure.');
    end
    signals=struct2cell(signal);
else
    signals={signal};
end
order=1:numel(domains);
if numel(domains)==2, order=[2 1]; end
for n=1:numel(signals)
    if ~isfloat(signals{n})||isempty(signals{n})
        error('signal components must be non-empty floating-point arrays.');
    end
    for k=1:numel(domains)
        if strcmp(domains{order(k)},'frequency')&&size(signals{n},k)<=2
            error('frequency dimensions require at least three points for ft_axis.');
        end
    end
end
end


