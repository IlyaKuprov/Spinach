% Builds JCAMP blocks from sampled magnetic resonance arrays. Syntax:
%
%         blocks=jcamp_grid(axes,signal,units,names,block)
%
% Parameters:
%
%    axes   - cell row of finite real floating-point column vectors,
%             in MATLAB dimension order; one to three dimensions
%
%    signal - floating-point array with size(signal,n)=numel(axes{n});
%             a 1D signal is a column. Alternatively, a scalar structure
%             of such arrays, e.g. pos/neg, cos/sin, or pos_pos/...;
%             all component arrays must have the declared shape
%
%    units  - cell row of JCAMP axis unit strings, one per dimension
%
%    names  - cell row of ASCII axis names, one per dimension
%
%    block  - scalar structure containing title, type, and metadata
%
% Outputs:
%
%    blocks - cell row of writer-ready blocks; named signal components
%             become separate linked blocks, identified by their names
%
% Arrays are sliced along dimension 1. Other dimensions become page
% coordinates. Real and imaginary values retain their signs; no FFT,
% quadrature recombination, conjugation, or normalisation is performed.
% This is the shared array-to-page translation used by the wrappers.
%
% talos@spindynamics.org

function blocks=jcamp_grid(axes,signal,units,names,block)

% Check consistency
grumble(axes,signal,units,names,block);

% Separate named quadrature or receiver components without recombination
if isstruct(signal)
    channels=fieldnames(signal);
    signals=struct2cell(signal);
else
    channels={''}; signals={signal};
end
blocks=cell(1,numel(signals));

% Translate each array into a trace or coordinate-labelled pages
for n=1:numel(signals)
    current=block; values=signals{n};
    if ~isempty(channels{n})
        current.title=[block.title '/' channels{n}];
        current.metadata=[current.metadata; {'$SPINACH COMPONENT',channels{n}}];
    end
    if isscalar(axes)
        current.x=axes{1}; current.y=values;
        current.xunits=units{1}; current.yunits='ARBITRARY UNITS';
        current.xname=names{1}; current.yname='Signal';
    else

        % Declare physical coordinates and ordinate components
        symbols=arrayfun(@(k) ['X' num2str(k)],1:numel(axes),'UniformOutput',false);
        variables=struct('name',names,'symbol',symbols,...
                         'type','INDEPENDENT','units',units);
        if isreal(values)
            components={values}; depend={'Y'}; labels={'Signal'};
        else
            components={real(values),imag(values)};
            depend={'R','I'}; labels={'Real signal','Imaginary signal'};
        end
        for k=1:numel(depend)
            variables(end+1)=struct('name',labels{k},'symbol',depend{k},...
                                   'type','DEPENDENT','units','ARBITRARY UNITS'); %#ok<AGROW>
        end
        current.variables=variables;

        % Traverse MATLAB columns in their original storage order
        lengths=cellfun(@numel,axes); npages=prod(lengths(2:end));
        pages=repmat(struct('x',[],'y',[],'xvar',symbols{1},...
                            'yvar','','coordinates',{{}}),1,npages*numel(depend));
        for k=1:npages
            coordinates=cell(numel(axes)-1,2); remainder=k-1;
            for m=2:numel(axes)
                index=mod(remainder,lengths(m))+1;
                remainder=floor(remainder/lengths(m));
                coordinates(m-1,:)={symbols{m},axes{m}(index)};
            end
            for m=1:numel(depend)
                page=(k-1)*numel(depend)+m;
                pages(page).x=axes{1};
                pages(page).y=components{m}((k-1)*lengths(1)+(1:lengths(1))).';
                pages(page).yvar=depend{m};
                pages(page).coordinates=coordinates;
            end
        end
        current.pages=pages;
    end
    blocks{n}=current;
end

end

% Consistency enforcement
function grumble(axes,signal,units,names,block)
if ~iscell(axes)||~isrow(axes)||isempty(axes)||numel(axes)>3
    error('axes must be a cell row containing one to three column vectors.');
end
if ~iscell(units)||~iscell(names)||~isrow(units)||~isrow(names)||...
   numel(units)~=numel(axes)||numel(names)~=numel(axes)
    error('units and names must be cell rows with one entry per axis.');
end
for n=1:numel(axes)
    if ~isfloat(axes{n})||~isreal(axes{n})||~iscolumn(axes{n})||...
       isempty(axes{n})||any(~isfinite(axes{n}))
        error('each axis must be a non-empty finite real floating-point column.');
    end
    if ~ischar(units{n})||~isrow(units{n})||isempty(units{n})||...
       ~ischar(names{n})||~isrow(names{n})||isempty(names{n})
        error('axis units and names must be non-empty character rows.');
    end
end
if ~isstruct(block)||~isscalar(block)||...
   ~all(isfield(block,{'title','type','metadata'}))
    error('block must contain title, type, and metadata.');
end
if isstruct(signal)
    if ~isscalar(signal)||isempty(fieldnames(signal))
        error('signal must be a non-empty scalar component structure.');
    end
    signals=struct2cell(signal);
else
    signals={signal};
end
for n=1:numel(signals)
    values=signals{n};
    if ~isfloat(values)||isempty(values)||any(isinf(values(:)))||...
       ndims(values)>max(2,numel(axes))
        error('signal components must be floating-point arrays without infinity.');
    end
    for k=1:numel(axes)
        if size(values,k)~=numel(axes{k})
            error('signal component dimensions must match the supplied axes.');
        end
    end
    if isscalar(axes)&&~iscolumn(values)
        error('a one-dimensional signal must be a column vector.');
    end
end
end


