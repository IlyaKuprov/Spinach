% Removes isotope numbers from nuclear isotope specifications. Syntax:
%
%                  elements=isot2elem(isotopes)
%
% Parameters:
%
%    isotopes  - cell array of character row vectors, such as
%                {'1H','13C','15N','35Cl'}
%
% Outputs:
%
%    elements  - cell array of element symbols with the same
%                shape and order as isotopes
%
% Only digits are removed; other characters are left unchanged.
% This is a string conversion, not an isotope-table lookup.
%
% talos@spindynamics.org

function elements=isot2elem(isotopes)

% Check consistency
grumble(isotopes);

% Strip isotope numbers
elements=regexprep(isotopes,'[0-9]','');

end

% Consistency enforcement
function grumble(isotopes)
if (~iscell(isotopes))||...
   (~all(cellfun(@(x)ischar(x)&&isrow(x),isotopes),'all'))
    error('isotopes must be a cell array of character row vectors.');
end
end


