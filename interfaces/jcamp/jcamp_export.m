% Exports NMR and EPR data as JCAMP-DX 5.01 text. Syntax:
%
%                       text=jcamp_export(data)
%
% Parameters:
%
%    data.title, data.origin, data.owner - non-empty ASCII row strings
%
%    data.blocks - non-empty cell vector of data block structures;
%                  each has title, type, and metadata (N-by-2 cell array
%                  of JCAMP label/value pairs, without ## or =)
%
%    data.filename - if present, also write the text to this filename
%
% Each block contains exactly one of the following representations:
%
%    x, y, xunits, yunits - equally sized column vectors and unit strings;
%       x is finite and real, y is real or complex and may contain NaN
%       (missing observations). Optional xname and yname label the axes.
%       Exact linspace axes with a finite first ordinate use XYDATA;
%       others use XYPOINTS. Complex y
%       uses NTUPLES, with separate real (R) and imaginary (I) pages.
%
%    variables, pages - general NTUPLES. variables is a structure vector
%       with name, symbol, type, and units fields. type is INDEPENDENT,
%       DEPENDENT, or PAGE. Symbols are unique uppercase identifiers.
%       pages is a structure vector with x, y (real column vectors),
%       xvar, yvar (variable symbols), and coordinates (N-by-2 cell array
%       of symbol/scalar pairs). Every variable occurs in the pages.
%       Each page specifies all remaining independent/page coordinates.
%       Pages may have different lengths and different x axes. This
%       represents arbitrary-dimensional, irregular, multichannel, and
%       hypercomplex data without guessing an array dimension order.
%
%    peaks, xunits, yunits - peaks is a scalar structure with x and
%       optionally y, width (real column vectors), multiplicity (cell
%       column of S,D,T,Q,M,U for NMR), assignment (cell column of ASCII strings),
%       and method (ASCII string). Unassigned widths/multiplicity require y;
%       EMR widths also require y. Assignments permit x-only NMR/EMR data.
%       method is required
%       for widths or EMR assignments, and describes the peak
%       finding/width convention. Width uses xunits. Assignment text
%       is enclosed in angle brackets; delimiters must not occur in it.
%
% Types: NMR FID, NMR SPECTRUM, NMR PEAK TABLE, NMR PEAK ASSIGNMENTS,
%        EMR SIMULATION, EMR MEASUREMENT. EPR uses the EMR data types.
% metadata holds experiment/sample information, including required
% technique-specific records; values are ASCII row strings, cell vectors
% of ASCII lines, or finite real numeric vectors. Examples:
% {'.OBSERVE FREQUENCY',400; '.OBSERVE NUCLEUS','^1H'} (MHz for NMR),
% {'.DETECTION MODE','CW'; '.METHOD','SPECTRUM'} for EMR. Nothing is
% inferred about nuclei, frequency, pulse sequence, units, or ownership.
% See README.md in this folder for metadata requirements and examples.
%
% Outputs:
%
%    text - ASCII character row containing the complete JCAMP file;
%           multiple blocks are enclosed in a LINK block
%
% Notes: floating-point AFFN has 17 significant digits and unit factors;
% no FFT, normalisation, quantisation, or unit conversion is performed.
% Infinity is rejected. NaN ordinates are written as ?; unavailable
% scalar statistics are omitted and NTUPLES attributes are left empty.
% Records are at most 80 characters. No third-party toolbox is needed. A requested
% file is replaced only after successful serialisation and writing.
%
% References: Davies and Lampen, Appl. Spectrosc. 47, 1093 (1993);
% Lampen et al., doi:10.1351/pac199971081549;
% Cammack et al., doi:10.1351/pac200678030613.
%
% talos@spindynamics.org

function text=jcamp_export(data)

% Check the complete input before opening a file
grumble(data);

% Create a LINK header for multiple independent datasets
lines={};
if numel(data.blocks)>1
    lines=[record('TITLE',data.title); record('JCAMP-DX','5.01');...
           record('DATA TYPE','LINK'); record('ORIGIN',data.origin);...
           record('OWNER',data.owner); record('BLOCKS',numel(data.blocks))];
end

% Serialise each independent dataset
for n=1:numel(data.blocks)
    block=data.blocks{n};
    lines=[lines; record('TITLE',block.title); record('JCAMP-DX','5.01');...
           record('DATA TYPE',block.type)]; %#ok<AGROW>

    % Select the representation from the supplied data
    if isfield(block,'peaks')
        data_class='PEAK TABLE';
        if isfield(block.peaks,'assignment'), data_class='PEAK ASSIGNMENTS'; end
    elseif isfield(block,'variables')||~isreal(block.y)
        data_class='NTUPLES';
    elseif numel(block.x)>1&&isfinite(block.y(1))&&isequal(block.x,linspace(block.x(1),block.x(end),numel(block.x))')
        data_class='XYDATA';
    else
        data_class='XYPOINTS';
    end
    class_label=data_class;
    if strcmp(data_class,'PEAK ASSIGNMENTS')&&startsWith(block.type,'NMR'), class_label='ASSIGNMENTS'; end
    lines=[lines; record('DATA CLASS',class_label); record('ORIGIN',data.origin);...
           record('OWNER',data.owner)]; %#ok<AGROW>
    if numel(data.blocks)>1, lines=[lines; record('BLOCK_ID',n)]; end %#ok<AGROW>

    % Preserve caller-supplied experiment and sample records
    for k=1:size(block.metadata,1)
        value=block.metadata{k,2};
        if startsWith(block.type,'NMR')&&strcmpi(regexprep(block.metadata{k,1},'[ _/-]',''),'.DELAY')&&isnumeric(value)
            value=['(' strjoin(arrayfun(@number,value(:)','UniformOutput',false),', ') ')'];
        end
        lines=[lines; record(block.metadata{k,1},value)]; %#ok<AGROW>
    end

    % Write general NTUPLES or split a complex trace into its two components
    if strcmp(data_class,'NTUPLES')
        if ~isfield(block,'variables')
            block.variables=struct('name',{'X','REAL','IMAGINARY','PAGE NUMBER'},...
                                   'symbol',{'X','R','I','N'},...
                                   'type',{'INDEPENDENT','DEPENDENT','DEPENDENT','PAGE'},...
                                   'units',{block.xunits,block.yunits,block.yunits,''});
            block.pages=struct('x',{block.x,block.x},'y',{real(block.y),imag(block.y)},...
                               'xvar',{'X','X'},'yvar',{'R','I'},...
                               'coordinates',{{'N',1},{'N',2}});
            if isfield(block,'xname'), block.variables(1).name=block.xname; end
            if isfield(block,'yname')
                block.variables(2).name=[block.yname '/REAL'];
                block.variables(3).name=[block.yname '/IMAGINARY'];
            end
        end
        lines=[lines; tuple_records(block)]; %#ok<AGROW>
    else

        % Write the common two-variable attributes
        if isfield(block,'peaks'), x=block.peaks.x; else, x=block.x; end
        lines=[lines; record('XUNITS',block.xunits); record('YUNITS',block.yunits);...
               record('XFACTOR',1); record('YFACTOR',1); record('NPOINTS',numel(x));...
               record('FIRSTX',x(1)); record('LASTX',x(end));...
               record('MINX',min(x)); record('MAXX',max(x))]; %#ok<AGROW>
        if isfield(block,'xname'), lines=[lines; record('XLABEL',block.xname)]; end %#ok<AGROW>
        if isfield(block,'yname'), lines=[lines; record('YLABEL',block.yname)]; end %#ok<AGROW>
        if isfield(block,'peaks')
            if isfield(block.peaks,'y'), y=block.peaks.y; else, y=[]; end
        else
            y=block.y;
        end
        if ~isempty(y)
            lines=[lines; record('FIRSTY',y(1)); record('MINY',min(y,[],'omitnan'));...
                   record('MAXY',max(y,[],'omitnan'))]; %#ok<AGROW>
        end

        % Write samples or a peak table without changing their order
        if ~isfield(block,'peaks')
            lines=[lines; table_records(x,y,'X','Y',data_class,data_class)]; %#ok<AGROW>
        else
            peaks=block.peaks; symbols='X';
            if isfield(peaks,'y'), symbols=[symbols 'Y']; end %#ok<AGROW>
            if isfield(peaks,'multiplicity'), symbols=[symbols 'M']; end %#ok<AGROW>
            if isfield(peaks,'width'), symbols=[symbols 'W']; end %#ok<AGROW>
            if isfield(peaks,'assignment'), symbols=[symbols 'A']; end %#ok<AGROW>
            if strcmp(data_class,'PEAK TABLE')&&startsWith(block.type,'NMR'), symbols=[symbols '..' symbols]; end %#ok<AGROW>
            lines=[lines; record(data_class,['(' symbols ')'])]; %#ok<AGROW>
            if isfield(peaks,'method'), lines=[lines; cellfun(@(s)['$$ ' s],wrap_record(peaks.method,77,0),'UniformOutput',false)]; end %#ok<AGROW>
            for k=1:numel(peaks.x)
                values={number(peaks.x(k))};
                if isfield(peaks,'y'), values{end+1}=number(peaks.y(k)); end %#ok<AGROW>
                if isfield(peaks,'multiplicity'), values{end+1}=peaks.multiplicity{k}; end %#ok<AGROW>
                if isfield(peaks,'width'), values{end+1}=number(peaks.width(k)); end %#ok<AGROW>
                if isfield(peaks,'assignment')
                    values{end+1}=['<' peaks.assignment{k} '>']; %#ok<AGROW>
                    row=['(' strjoin(values,', ') ')'];
                else
                    row=strjoin(values,', ');
                end
                lines=[lines; wrap_record(row,80,0)]; %#ok<AGROW>
            end
        end
    end
    lines=[lines; record('END','')]; %#ok<AGROW>
end
if numel(data.blocks)>1, lines=[lines; record('END','')]; end
text=[strjoin(lines',char([13 10])) char([13 10])];

% Replace the destination only after the complete file has been written
if isfield(data,'filename')
    folder=fileparts(data.filename);
    if isempty(folder), folder=pwd; end
    temporary=tempname(folder); file_id=fopen(temporary,'w');
    if file_id<0, error('cannot open a temporary JCAMP file.'); end
    try
        written=fwrite(file_id,text,'char');
        status=fclose(file_id); file_id=-1;
        if written~=numel(text)||status~=0, error('could not write the complete JCAMP file.'); end
        [success,message]=movefile(temporary,data.filename,'f');
        if ~success, error('could not replace the JCAMP file: %s',message); end
    catch exception
        if file_id>=0, fclose(file_id); end
        if isfile(temporary), delete(temporary); end
        rethrow(exception);
    end
end

end

% Serialise shared variable attributes and individually sampled NTUPLES pages
function lines=tuple_records(block)

% Derive attributes from each variable's actual occurrences
variables=block.variables; pages=block.pages;
count=numel(variables); dimensions=zeros(1,count);
first=zeros(1,count); last=zeros(1,count); low=zeros(1,count); high=zeros(1,count);
for n=1:count
    values={}; coords=[]; dimension=0;
    for k=1:numel(pages)
        if strcmp(variables(n).symbol,pages(k).xvar)
            values{end+1}=double(pages(k).x); dimension=max(dimension,numel(pages(k).x)); %#ok<AGROW>
        elseif strcmp(variables(n).symbol,pages(k).yvar)
            values{end+1}=double(pages(k).y); dimension=max(dimension,numel(pages(k).y)); %#ok<AGROW>
        else
            idx=strcmp(variables(n).symbol,pages(k).coordinates(:,1));
            if any(idx)
                coords(end+1)=pages(k).coordinates{idx,2}; %#ok<AGROW>
                values{end+1}=double(pages(k).coordinates{idx,2}); %#ok<AGROW>
            end
        end
    end
    sample=vertcat(values{:});
    if dimension==0, dimensions(n)=numel(unique(coords)); else, dimensions(n)=dimension; end
    first(n)=sample(1); last(n)=sample(end);
    low(n)=min(sample,[],'omitnan'); high(n)=max(sample,[],'omitnan');
end

% Write the NTUPLES attribute table with lossless AFFN scale factors
lines=[record('NTUPLES',block.type); record('VAR_NAME',strjoin({variables.name},', '));...
       record('SYMBOL',strjoin({variables.symbol},', '));...
       record('VAR_TYPE',strjoin({variables.type},', '));...
       record('VAR_FORM',strjoin(repmat({'AFFN'},1,count),', '));...
       record('VAR_DIM',dimensions); record('UNITS',strjoin({variables.units},', '));...
       record('FIRST',first); record('LAST',last); record('MIN',low);...
       record('MAX',high); record('FACTOR',ones(1,count))];

% Preserve page coordinates, point counts, and unequal sampling
for n=1:numel(pages)
    page=pages(n); coords=cell(size(page.coordinates,1),1);
    for k=1:numel(coords)
        coords{k}=[page.coordinates{k,1} '=' number(page.coordinates{k,2})];
    end
    lines=[lines; record('PAGE',strjoin(coords,', ')); record('NPOINTS',numel(page.x))]; %#ok<AGROW>
    xidx=find(strcmp(page.xvar,{variables.symbol})); yidx=strcmp(page.yvar,{variables.symbol});
    if numel(page.x)>1&&isfinite(page.y(1))&&isfinite(first(yidx))&&numel(page.x)==dimensions(xidx)&&...
       page.x(1)==first(xidx)&&page.x(end)==last(xidx)&&...
       isequal(page.x,linspace(page.x(1),page.x(end),numel(page.x))')
        data_class='XYDATA';
    else
        data_class='XYPOINTS';
    end
    lines=[lines; table_records(page.x,page.y,page.xvar,page.yvar,data_class,'DATA TABLE')]; %#ok<AGROW>
end
lines=[lines; record('END NTUPLES',block.type)];

end

% Serialise regular or explicitly tabulated pairs with no numerical rescaling
function lines=table_records(x,y,xvar,yvar,data_class,label)

% Select an incremental or an explicit-pair variable list
if strcmp(data_class,'XYDATA')
    descriptor=['(' xvar '++(' yvar '..' yvar '))'];
else
    descriptor=['(' xvar yvar '..' xvar yvar ')'];
end
if strcmp(label,'DATA TABLE')
    if strcmp(data_class,'XYDATA'), descriptor=[descriptor ', XYDATA']; else, descriptor=[descriptor ', PROFILE']; end
end
lines=record(label,descriptor);

% One pair per line preserves every explicit abscissa and fits the record limit
rows=cell(numel(x),1);
for n=1:numel(x)
    if strcmp(data_class,'XYDATA')
        rows{n}=[number(x(n)) ' ' number(y(n))];
    else
        rows{n}=[number(x(n)) ', ' number(y(n))];
    end
end
lines=[lines; rows];

end

% Format a numeric scalar with sufficient digits to round-trip an IEEE double
function text=number(value)

% Use the JCAMP missing-value marker instead of a non-standard NaN token
if isnan(value)
    text='?';
elseif isinteger(value)
    if startsWith(class(value),'uint'), text=sprintf('%u',value); else, text=sprintf('%d',value); end
else
    text=sprintf('%.17G',value);
end

end

% Serialise an ASCII labelled data record and preserve multiline values
function lines=record(label,value)

% Convert numeric vectors into comma-separated AFFN values
if isnumeric(value)
    if isscalar(value)&&isnan(value), lines=cell(0,1); return; end
    values=arrayfun(@number,value(:)','UniformOutput',false);
    values(isnan(value(:)'))={''};
    value=strjoin(values,', ');
end
if ischar(value), value={value}; end
prefix=['##' upper(label) '='];
lines=wrap_record([prefix value{1}],80,numel(prefix));
for n=2:numel(value), lines=[lines; wrap_record(value{n},80,0)]; end %#ok<AGROW>

end

% Wrap on existing spaces without splitting numbers, symbols, or assignment tags
function lines=wrap_record(value,limit,protected)

% Continuation records have no label prefix
lines={};
while numel(value)>limit
    split=find(value(1:limit)==' '&(1:limit)>protected,1,'last');
    if isempty(split)&&protected>0
        lines{end+1,1}=value(1:protected); value=value(protected+1:end); protected=0; %#ok<AGROW>
        continue;
    end
    if isempty(split), error('a JCAMP token exceeds the available record width.'); end
    lines{end+1,1}=value(1:split-1); value=value(split+1:end); protected=0; %#ok<AGROW>
end
lines{end+1,1}=value;

end

% Recognise ASCII text that cannot introduce a labelled record or a comment
function valid=ascii_text(value)

% Empty unit strings are permitted for page counters
valid=ischar(value)&&(isrow(value)||isempty(value))&&...
      all(value>=32&value<=126)&&~contains(value,'##')&&~contains(value,'$$');

end

% Recognise a non-empty floating-point column without converting its shape
function valid=numeric_column(value)

% Floating-point inputs avoid silent conversion of large integer samples
valid=isfloat(value)&&iscolumn(value)&&~isempty(value)&&~issparse(value);

end

% Check caller-controlled representations, units, metadata, and page coordinates
function grumble(data)
if ~isstruct(data)||~isscalar(data), error('data must be a scalar structure.'); end
for field={'title','origin','owner'}
    if ~isfield(data,field{1})||~ascii_text(data.(field{1}))||isempty(strtrim(data.(field{1})))
        error('data.%s must be a non-empty ASCII row string.',field{1});
    end
end
if isfield(data,'filename')&&(~ischar(data.filename)||~isrow(data.filename)||isempty(data.filename)||isfolder(data.filename))
    error('data.filename must be a non-empty character row naming a file, not a directory.');
end
if ~isfield(data,'blocks')||~iscell(data.blocks)||~isvector(data.blocks)||isempty(data.blocks)
    error('data.blocks must be a non-empty cell vector.');
end
reserved={'TITLE','JCAMPDX','DATATYPE','DATACLASS','ORIGIN','OWNER','BLOCKS','BLOCKID',...
          'END','XUNITS','YUNITS','XFACTOR','YFACTOR','FIRSTX','LASTX','FIRSTY',...
          'MINX','MAXX','MINY','MAXY','NPOINTS','DELTAX','XLABEL','YLABEL','XYDATA',...
          'XYPOINTS','PEAKTABLE','PEAKASSIGNMENTS','NTUPLES','ENDNTUPLES','VARNAME',...
          'SYMBOL','VARTYPE','VARFORM','VARDIM','UNITS','FIRST','LAST','MIN','MAX',...
          'FACTOR','PAGE','DATATABLE'};
for n=1:numel(data.blocks)
    block=data.blocks{n};
    if ~isstruct(block)||~isscalar(block)||~isfield(block,'title')||...
       ~ascii_text(block.title)||isempty(strtrim(block.title))
        error('each block needs a non-empty ASCII title.');
    end
    types={'NMR FID','NMR SPECTRUM','NMR PEAK TABLE','NMR PEAK ASSIGNMENTS',...
           'EMR SIMULATION','EMR MEASUREMENT'};
    if ~isfield(block,'type')||~ischar(block.type)||~ismember(block.type,types)
        error('block.type must be a documented NMR or EMR data type.');
    end
    if ~isfield(block,'metadata')||~iscell(block.metadata)||size(block.metadata,2)~=2||~ismatrix(block.metadata)
        error('block.metadata must be an N-by-2 cell array.');
    end
    labels=cell(size(block.metadata,1),1);
    for k=1:size(block.metadata,1)
        label=block.metadata{k,1}; value=block.metadata{k,2};
        if ~ischar(label)||isempty(regexp(label,'^[.$]?[A-Za-z][A-Za-z0-9 _/.-]*$','once'))||numel(label)>77
            error('metadata labels must be ASCII JCAMP identifiers, without ## or =.');
        end
        if ~startsWith(label,'$')&&contains(label(2:end),'.')
            error('periods may only prefix technique-specific labels; use $ for private labels.');
        end
        labels{k}=upper(regexprep(label,'[ _/-]',''));
        if ismember(labels{k},reserved), error('metadata must not override generated label %s.',label); end
        if ischar(value)
            valid=ascii_text(value);
        elseif iscell(value)
            valid=isvector(value)&&~isempty(value)&&all(cellfun(@ascii_text,value));
        else
            valid=isnumeric(value)&&isreal(value)&&isvector(value)&&~isempty(value)&&all(isfinite(value));
        end
        if ~valid, error('metadata values must be ASCII text or finite real numeric vectors.'); end
    end
    if numel(unique(labels))~=numel(labels), error('metadata labels must be unique after JCAMP canonicalisation.'); end
    if startsWith(block.type,'NMR')
        required={'.OBSERVEFREQUENCY','.OBSERVENUCLEUS'};
        if strcmp(block.type,'NMR FID'), required=[required {'.DELAY','.ACQUISITIONMODE'}]; end %#ok<AGROW>
    else
        required={'.DETECTIONMODE','.METHOD'};
        if strcmp(block.type,'EMR SIMULATION')
            required=[required {'.SIMULATIONSOURCE','.SIMULATIONPARAMETERS'}]; %#ok<AGROW>
        end
    end
    if ~all(ismember(required,labels)), error('missing required technique metadata: %s.',strjoin(required,', ')); end
    for k=find(ismember(labels,required))'
        value=block.metadata{k,2};
        if isempty(value)||(ischar(value)&&isempty(strtrim(value)))||...
           (iscell(value)&&all(cellfun(@(s)isempty(strtrim(s)),value)))
            error('required technique metadata must not be empty.');
        end
    end
    if startsWith(block.type,'NMR')
        value=block.metadata{strcmp(labels,'.OBSERVEFREQUENCY'),2};
        if ischar(value), value=str2double(value); end
        if ~isnumeric(value)||~isscalar(value)||~isreal(value)||~isfinite(value)||value<=0
            error('.OBSERVE FREQUENCY must be a positive frequency in MHz.');
        end
        value=block.metadata{strcmp(labels,'.OBSERVENUCLEUS'),2};
        if ~ischar(value)||isempty(regexp(value,'^\^[0-9]+[A-Z][a-z]?$','once'))
            error('.OBSERVE NUCLEUS must use a JCAMP isotope label such as ^1H.');
        end
        idx=strcmp(labels,'.DELAY');
        if any(idx)
            value=block.metadata{idx,2};
            if ischar(value)
                values=regexp(value,'^\(\s*([^,]+),\s*([^,]+)\)$','tokens','once');
                values=str2double(values);
                valid=numel(values)==2&&isreal(values)&&all(isfinite(values));
            else
                valid=isnumeric(value)&&numel(value)==2;
            end
            if ~valid, error('.DELAY must contain two finite real pre-acquisition delays, numeric or (RD, ID) text.'); end
        end
    end
    if startsWith(block.type,'EMR')
        methods={'DYNAMIC','ELDOR','ENDOR','ESEEM','ODMR','GONIOMETER','HYSCORE',...
                 'KINETIC','SATURATION','SPECTRUM','FID','TRIPLE','IMAGING','SPECTRAL SPATIAL'};
        value=block.metadata{strcmp(labels,'.METHOD'),2};
        if ~ischar(value)||~ismember(value,methods), error('.METHOD must be a documented EMR method identifier.'); end
        for k=find(ismember(labels,{'.SIMULATIONSOURCE','.SIMULATIONPARAMETERS'}))'
            if ~ischar(block.metadata{k,2})&&~iscell(block.metadata{k,2})
                error('EMR simulation source and parameters must contain ASCII text.');
            end
        end
    end
    modes={'.ACQUISITIONMODE',{'SIMULTANEOUS','SEQUENTIAL','SINGLE'}; '.DETECTIONMODE',{'CW','PULSE'}};
    for k=1:size(modes,1)
        idx=strcmp(modes{k,1},labels);
        if any(idx)&&(~ischar(block.metadata{idx,2})||~ismember(block.metadata{idx,2},modes{k,2}))
            error('invalid %s metadata.',modes{k,1});
        end
    end
    choices=[isfield(block,'x')||isfield(block,'y'),...
             isfield(block,'variables')||isfield(block,'pages'),isfield(block,'peaks')];
    if sum(choices)~=1, error('each block must supply exactly one of x/y, variables/pages, or peaks.'); end
    if strcmp(block.type,'NMR FID')
        allowed_units={'SECONDS'};
    elseif strcmp(block.type,'NMR SPECTRUM')
        allowed_units={'HZ'};
    elseif startsWith(block.type,'NMR')
        allowed_units={'HZ','PPM'};
    else
        allowed_units={'DEGREE','HERTZ','KELVIN','SECOND','TESLA','WATT'};
    end
    if startsWith(block.type,'NMR')
        allowed_yunits={'ARBITRARY UNITS','MAGNITUDE','POWER'};
    else
        allowed_yunits={'ARBITRARY UNITS','INTENSITY','POWER'};
    end
    if ~choices(2)
        for field={'xunits','yunits'}
            if ~isfield(block,field{1})||~ascii_text(block.(field{1}))||isempty(strtrim(block.(field{1})))||contains(block.(field{1}),',')
                error('each two-variable block needs non-empty ASCII xunits and yunits without commas.');
            end
        end
        if ~ismember(block.xunits,allowed_units), error('xunits must match the declared NMR or EMR data type.'); end
        if ~ismember(block.yunits,allowed_yunits), error('yunits must match the declared NMR or EMR data type.'); end
        for field={'xname','yname'}
            if isfield(block,field{1})&&(~ascii_text(block.(field{1}))||isempty(block.(field{1}))||contains(block.(field{1}),','))
                error('axis names must be non-empty ASCII row strings without commas.');
            end
        end
    end
    if choices(1)
        if ~isfield(block,'x')||~isfield(block,'y')||~numeric_column(block.x)||...
           ~isreal(block.x)||any(~isfinite(block.x))||~numeric_column(block.y)||...
           numel(block.x)~=numel(block.y)||any(isinf(real(block.y))|isinf(imag(block.y)))
            error('x/y must be equally sized floating-point columns; x finite real, y finite or missing.');
        end
        if startsWith(block.type,'NMR')&&~isreal(block.y)&&~strcmp(block.yunits,'ARBITRARY UNITS')
            error('complex NMR traces require ARBITRARY UNITS; supply real transformed magnitude or power data.');
        end
        if contains(block.type,'PEAK'), error('NMR peak data types require a peaks structure.'); end
    elseif choices(2)
        if contains(block.type,'PEAK'), error('NMR peak data types require a peaks structure.'); end
        if ~isfield(block,'variables')||~isstruct(block.variables)||~isvector(block.variables)||isempty(block.variables)||...
           ~all(isfield(block.variables,{'name','symbol','type','units'}))
            error('variables must be a non-empty structure vector with name, symbol, type, and units.');
        end
        variables=block.variables; symbols={variables.symbol};
        for k=1:numel(variables)
            if ~ascii_text(variables(k).name)||isempty(variables(k).name)||contains(variables(k).name,',')||...
               ~ascii_text(variables(k).units)||contains(variables(k).units,',')||...
               ~ischar(symbols{k})||isempty(regexp(symbols{k},'^[A-Z][A-Z0-9]*$','once'))||...
               ~ischar(variables(k).type)||~ismember(variables(k).type,{'INDEPENDENT','DEPENDENT','PAGE'})
                error('invalid NTUPLES variable descriptor.');
            end
        end
        for k=1:numel(variables)
            if ~strcmp(variables(k).type,'PAGE')&&isempty(strtrim(variables(k).units))
                error('independent and dependent variables require explicit units.');
            end
        end
        if numel(unique(symbols))~=numel(symbols), error('variable symbols must be unique.'); end
        if ~isfield(block,'pages')||~isstruct(block.pages)||~isvector(block.pages)||isempty(block.pages)||...
           ~all(isfield(block.pages,{'x','y','xvar','yvar','coordinates'}))
            error('pages must supply x, y, xvar, yvar, and coordinates.');
        end
        seen=false(size(symbols));
        for k=1:numel(block.pages)
            page=block.pages(k);
            if ~numeric_column(page.x)||~isreal(page.x)||any(~isfinite(page.x))||...
               ~numeric_column(page.y)||~isreal(page.y)||any(isinf(page.y))||numel(page.x)~=numel(page.y)
                error('page x/y must be equally sized real floating-point columns; x finite, y finite or missing.');
            end
            xidx=find(strcmp(page.xvar,symbols)); yidx=find(strcmp(page.yvar,symbols));
            if numel(xidx)~=1||numel(yidx)~=1||~strcmp(variables(xidx).type,'INDEPENDENT')||...
               ~strcmp(variables(yidx).type,'DEPENDENT')
                error('page xvar/yvar must name an independent/dependent variable pair.');
            end
            if ~ismember(variables(xidx).units,allowed_units)
                error('tabulated abscissa units must match the declared NMR or EMR data type.');
            end
            if ~ismember(variables(yidx).units,allowed_yunits)
                error('tabulated ordinate units must match the declared NMR or EMR data type.');
            end
            if ~iscell(page.coordinates)||size(page.coordinates,2)~=2||~ismatrix(page.coordinates)
                error('page.coordinates must be an N-by-2 cell array.');
            end
            fixed=false(size(symbols));
            for m=1:size(page.coordinates,1)
                idx=find(strcmp(page.coordinates{m,1},symbols)); value=page.coordinates{m,2};
                if numel(idx)~=1||idx==xidx||strcmp(variables(idx).type,'DEPENDENT')||fixed(idx)||...
                   ~isnumeric(value)||~isscalar(value)||~isreal(value)||~isfinite(value)||...
                   (isinteger(value)&&abs(value)>flintmax)
                    error('page coordinates must uniquely fix non-tabulated variables with double-representable scalars.');
                end
                fixed(idx)=true;
            end
            expected=~strcmp({variables.type},'DEPENDENT'); expected(xidx)=false;
            if ~isequal(fixed,expected), error('each page must fix all remaining independent/page variables.'); end
            seen=seen|fixed; seen([xidx yidx])=true;
        end
        if ~all(seen), error('every NTUPLES variable must occur in the pages.'); end
    else
        peaks=block.peaks;
        if ~isstruct(peaks)||~isscalar(peaks)||~isfield(peaks,'x')||...
           ~numeric_column(peaks.x)||~isreal(peaks.x)||any(~isfinite(peaks.x))
            error('peaks.x must be a finite real floating-point column.');
        end
        for field={'y','width'}
            if isfield(peaks,field{1})&&(~numeric_column(peaks.(field{1}))||~isreal(peaks.(field{1}))||...
               any(~isfinite(peaks.(field{1})))||numel(peaks.(field{1}))~=numel(peaks.x))
                error('peak numeric columns must be finite real and match peaks.x.');
            end
        end
        if isfield(peaks,'width')&&(any(peaks.width<0)||(~isfield(peaks,'y')&&...
           (~isfield(peaks,'assignment')||startsWith(block.type,'EMR'))))
            error('widths must be non-negative; unassigned and EMR peaks require y.');
        end
        for field={'multiplicity','assignment'}
            if isfield(peaks,field{1})&&(~iscell(peaks.(field{1}))||~iscolumn(peaks.(field{1}))||...
               numel(peaks.(field{1}))~=numel(peaks.x)||~all(cellfun(@ascii_text,peaks.(field{1}))))
                error('peak text columns must be ASCII cell columns matching peaks.x.');
            end
        end
        if isfield(peaks,'multiplicity')&&((~isfield(peaks,'y')&&~isfield(peaks,'assignment'))||...
           ~all(ismember(peaks.multiplicity,{'S','D','T','Q','M','U'})))
            error('unassigned multiplicity requires y; NMR multiplicity uses S,D,T,Q,M,U.');
        end
        if startsWith(block.type,'EMR')&&isfield(peaks,'multiplicity'), error('multiplicity is defined only for NMR peaks.'); end
        if ~isfield(peaks,'assignment')&&isfield(peaks,'width')&&isfield(peaks,'multiplicity')
            error('unassigned peak tables support either width or multiplicity, not both.');
        end
        if isfield(peaks,'assignment')&&any(cellfun(@(s)~isempty(regexp(s,'[<>\r\n]','once')),peaks.assignment))
            error('assignment strings must not contain angle brackets or newlines.');
        end
        if (isfield(peaks,'width')||(isfield(peaks,'assignment')&&startsWith(block.type,'EMR')))&&~isfield(peaks,'method')
            error('peak widths and EMR assignments require a method description.');
        end
        if isfield(peaks,'method')&&(~ascii_text(peaks.method)||isempty(strtrim(peaks.method)))
            error('peak method must be a non-empty ASCII row string.');
        end
        if ~isfield(peaks,'assignment')&&~isfield(peaks,'y'), error('unassigned peak tables require y.'); end
        if startsWith(block.type,'NMR')&&(~contains(block.type,'PEAK')||...
           (isfield(peaks,'assignment')~=strcmp(block.type,'NMR PEAK ASSIGNMENTS')))
            error('NMR peak data type must agree with the presence of assignments.');
        end
    end
end
end


