% Export a complete CCE result without discarding complex phase. Syntax:
%
%                      cce_export(result,filename,format)
%
% format is 'mat','json', or 'csv'. MAT/JSON store a versioned envelope;
% JSON preserves numeric type, shape, complex components and nonfinite
% values. CSV stores numeric vectors/matrices as long-form real/imag rows.
% A JSON sidecar, filename.json, retains the full result and diagnostics.
% CSV is a human-readable view; its sidecar is the round-trip authority.
%
% ilya.kuprov@weizmann.ac.il

function cce_export(result,filename,format)

% Validate serialization choices before opening a destination
 grumble(result,filename,format);
 envelope=struct('schema','spinach-cce','version',1,'payload',result);

% Preserve native MATLAB arrays or explicitly typed JSON values
switch format
    case 'mat'
        save(filename,'envelope');
    case 'json'
        fid=fopen(filename,'w');
        if fid<0, error('Cannot open export destination.'); end
        cleanup=onCleanup(@()fclose(fid));
        fwrite(fid,jsonencode(pack(envelope)),'char');
    case 'csv'
        fid=fopen(filename,'w');
        if fid<0, error('Cannot open export destination.'); end
        cleanup=onCleanup(@()fclose(fid));
        fprintf(fid,'field,index,real,imag\n');
        names=fieldnames(result);
        for k=1:numel(names)
            value=result.(names{k});
            if isnumeric(value)
                for j=1:numel(value)
                    fprintf(fid,'%s,%d,%.17g,%.17g\n',names{k},j,real(value(j)),imag(value(j)));
                end
            end
        end
        cce_export(result,[filename '.json'],'json');
end

end

function encoded=pack(value)

% Encode type and dimensions instead of relying on lossy JSON inference
if isinteger(value)
    encoded=struct('kind','integer','class',class(value),'shape',size(value),...
                   'bytes',typecast(value(:),'uint8').');
elseif isnumeric(value)||islogical(value)
    encoded=struct('kind','array','class',class(value),'shape',size(value),...
                   'real',{arrayfun(@(v)sprintf('%.17g',v),real(value(:)),'UniformOutput',false)},...
                   'imag',{arrayfun(@(v)sprintf('%.17g',v),imag(value(:)),'UniformOutput',false)});
elseif ischar(value)
    encoded=struct('kind','char','value',value);
elseif isstring(value)
    encoded=struct('kind','string','shape',size(value),'items',{cellstr(value(:))});
elseif iscell(value)
    encoded=struct('kind','cell','shape',size(value),'items',{cellfun(@pack,value(:),'UniformOutput',false)});
elseif isstruct(value)
    names=fieldnames(value); entries=cell(numel(value),numel(names));
    for k=1:numel(value)
        for j=1:numel(names), entries{k,j}=pack(value(k).(names{j})); end
    end
    encoded=struct('kind','struct','shape',size(value),'fields',{names},'items',{entries(:)});
else
    error('Unsupported export type: %s.',class(value));
end

end

function grumble(result,filename,format)

% Require a scalar result structure and a supported file representation
if ~isstruct(result)||~isscalar(result)||~ischar(filename)||isempty(filename)||...
   ~ischar(format)||~ismember(format,{'mat','json','csv'})
    error('Invalid CCE export arguments.');
end

end
