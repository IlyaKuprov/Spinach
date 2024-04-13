% Reduces streak artefacts in 2D and 3D NMR spectra. Edges of the
% input spectrum must be free of genuine signals. Cell arrays and
% structures are processed recursively. Syntax:
%
%                  spectrum=destreak(spectrum)
%
% The function works by subtracting the kronecker propduct of edge
% lines from the spectrum matrix.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=Destreak.m>

function spectrum=destreak(spectrum)

% Process structures and cell arrays recursively
if isstruct(spectrum)
    
    % Get the field names
    struct_fieldnames=fieldnames(spectrum);
    
    % Loop over field names
    for n=1:length(struct_fieldnames)
        
        % Call itself for each field name
        spectrum.(struct_fieldnames{n})=destreak(spectrum.(struct_fieldnames{n}));
                      
    end
    
elseif iscell(spectrum)
    
    % Loop over cells
    parfor n=1:numel(spectrum)
        
        % Call itself for each cell
        spectrum{n}=destreak(spectrum{n});
        
    end
    
end

% Check consistency
grumble(spectrum);

% Decide problem dimensionality
switch ndims(spectrum)
    
    case 2
        
        % Destreak the spectrum
        spectrum=spectrum-kron(spectrum(:,1),ones(1,size(spectrum,2)));
        spectrum=spectrum-kron(spectrum(1,:),ones(size(spectrum,1),1));
        
    case 3
        
        % Destreak the spectrum
        spectrum=spectrum-repmat(spectrum(:,1,:),1,size(spectrum,2),1);
        spectrum=spectrum-repmat(spectrum(1,:,:),size(spectrum,1),1,1);
        spectrum=spectrum-repmat(spectrum(:,:,1),1,1,size(spectrum,3));
        
    otherwise
        
        % Complain and bomb out
        error('unsupported spectrum dimensionality.');
        
end

end

% Consistency enforcement
function grumble(spectrum)
if ~isnumeric(spectrum)
    error(['spectrum must be a numeric array, a cell array '... 
           'of numeric arrays or a structure with numeric arrays inside.']);
end
if isvector(spectrum)
    error('destreaking is not applicable to one-dimensional spectra.');
end
end

% If it flies, floats or fucks, you are better off renting it.
%
% Felix Dennis

