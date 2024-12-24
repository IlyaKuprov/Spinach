% Zeroes out the diagonal of a 2D spectrum using the brush 
% with the specified dimensions. Syntax:
%
%               spec=killdiag(spec,brush_dim)
%
% Parameters:
%
%     spec      - 2D matrix representing a spectrum
%
%     brush_dim - the width of the band to zero out
%                 around the diagonal, points
%
% Outputs:
%
%     spec      - 2D matrix representing a spectrum
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=killdiag.m>

function spec=killdiag(spec,brush_dim)

% Check consistency
grumble(spec,brush_dim);

% Loop over the column index
for n=1:size(spec,2)
    
    % Find the row index
    k=n*size(spec,1)/size(spec,2);
    
    % Find the row index extents
    k=floor(k-brush_dim/2):ceil(k+brush_dim/2);
    
    % Avoid array boundaries
    k(k<1)=[]; k(k>size(spec,1))=[];
    
    % Zero the elements
    spec(k,n)=0;
    
end

end

% Consistency enforcement
function grumble(spec,brush_dim)
if (~isnumeric(spec))||(~ismatrix(spec))
    error('spec must be a matrix.');
end
if (~isnumeric(brush_dim))||(~isscalar(brush_dim))||...
   (~isreal(brush_dim))||(brush_dim<1)||(mod(brush_dim,1)~=0)
    error('brush_dim must be a positive real integer.');
end
if any(size(spec)<brush_dim)
    error('the brush cannot be wider than the spectrum.');
end
end

% "I do wish we could chat longer, but... I'm having
%  an old friend for dinner."
%
% Hannibal Lecter

