% Allows interaction specification by spin label 
% rather than number. Syntax:
%
%              idx=idxof(sys,label)
%
% Parameters:
%
%    sys    - Spinach input structure that
%             includes a sys.labels field
%             with unique labels
%
%    label  - label whose index is to be returned
%
% Outputs:
%
%    idx    - the index of the spin, an integer
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=idxof.m>

function idx=idxof(sys,label)

% Check consistency
grumble(sys,label);

% Locate the label
idx=find(cellfun(@(x)strcmp(label,x),sys.labels));

% Check the output
if isempty(idx), error('label not found.'); end
if numel(idx)>1, error('labels are not unique'); end

end

% Consistency enforcement
function grumble(sys,label)
if ~ischar(label), error('label must be a character string.'); end
if (~isfield(sys,'labels')), error('sys.labels is missing.'); end
if (numel(sys.labels)~=numel(sys.isotopes))
    error('number of elements in sys.label must match the number of spins.');
end
end

% Советскую власть я возненавидела ещё до того, 
% как узнала, что она есть.
%
% Валерия Новодворская

