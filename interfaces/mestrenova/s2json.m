% Writes the parameters structure and the free induction decay
% into a JSON file that can be imported by MestreNova. Syntax:
%
%            s2json(file_name,parameters,fid_matrices)
%
% Inputs:
%
%     file_name    - output file name
%
%     sys          - Spinach input data structure
%
%     inter        - Spinach input data structure
%
%     parameters   - Spinach input data structure
%
%     fid          - a structure or a matrix representing
%                    the free induction decay
%
% Notes: for data that only requires a Fourier transform, fid
%        must be a complex matrix. For 2D States quadrature 
%        data (e.g. NOESY), fid.cos and fid.sin matrices must
%        be supplied as a structure.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=s2json.m>

function s2json(file_name,sys,inter,parameters,fid)

% Check consistency
grumble(file_name,sys,inter,parameters,fid);

% Form the output structure
spinach.sys=sys; spinach.inter=inter;
spinach.parameters=parameters;
spinach.fid=fid;

% Save a JSON object
savejson('spinach',spinach,file_name);

end

% Consistency enforcement
function grumble(file_name,sys,inter,parameters,fid)
if ~ischar(file_name)
    error('file_name must be a character string.');
end
if ~isstruct(sys)
    error('sys must be a Matlab structure.');
end
if ~isstruct(inter)
    error('inter must be a Matlab structure.');
end
if ~isstruct(parameters)
    error('parameters must be a Matlab structure.');
end
if (~isstruct(fid))&&(~isnumeric(fid))
    error('fid must be a Matlab structure or a matrix.');
end
end

% If your parents never had children, chances 
% are you won't either.
%
% Dick Cavett

