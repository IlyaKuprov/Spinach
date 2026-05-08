% transforms Liouville ket to Hilbert space matrix
%
% input parameters:
% Ml: the Liouville space ket
%
% output parameters:
% Mhilb: the Hilbert space matrix
%
% February 2024 A. Kehl (akehl@gwdg.de)
% Spinach functional refactor May 2026 Talos


function Mhilb=kehl_lket_to_mat(Ml)

    % Check consistency
    grumble(Ml);
hilb_size=sqrt(size(Ml,1));
Mhilb=reshape(Ml,hilb_size,hilb_size);

end

function grumble(Ml)
if ~isnumeric(Ml)
    error('Ml must be numeric.');
end
end

