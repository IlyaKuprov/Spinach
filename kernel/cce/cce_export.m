% Export a complete CCE result without discarding complex phase. Syntax:
%
%                      cce_export(filename,model,result)
%
% filename must end in .mat, .json, or .csv; model and result are scalar
% structures. MAT/JSON store a versioned envelope with both structures;
% JSON preserves numeric type, shape, complex components and nonfinite
% values. CSV stores numeric vectors/matrices as long-form real/imag rows.
% A JSON sidecar, filename.json, retains the full result and diagnostics.
% CSV is a human-readable view; its sidecar is the round-trip authority.
%
% Inputs:
%   filename - .mat, .json or .csv destination.
%   model, result - scalar structures to serialize.
%
% Outputs:
%   No return value; writes the requested file (and JSON sidecar for CSV).
%
% See also cce_import.
%
% ilya.kuprov@weizmann.ac.il

function cce_export(filename,model,result)

% Validate serialization choices before opening a destination
grumble(filename,model,result);

% Infer the representation without extending the fixed public signature
[~,~,extension]=fileparts(filename); format=lower(extension(2:end));
envelope=struct('schema','spinach-cce','version',1,'model',model,'payload',result);

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
        cce_export([filename '.json'],model,result);
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

function grumble(filename,model,result)

if ~ischar(filename)||~isrow(filename)||isempty(filename)||...
   ~isstruct(model)||~isscalar(model)||~isstruct(result)||~isscalar(result)
    error('Invalid CCE export arguments.');
end

[~,~,extension]=fileparts(filename);
if ~ismember(lower(extension),{'.mat','.json','.csv'})
    error('CCE export requires a .mat, .json or .csv filename.');
end

end


