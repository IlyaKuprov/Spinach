% transforms Hilbert space matrix to Liouville Bra using Spinach
%
% input parameters:
% Mh: the Hilbert space matrix
%
% output parameters:
% Ml: the Liouville space Ket
%
% February 2024 A. Kehl (akehl@gwdg.de)
% Spinach migration May 2026 Talos


function Ml=kehl_mat_to_lbra(Mhilb)

    % Check consistency
    grumble(Mhilb);
    Ml=full(hilb2liouv(sparse(Mhilb.'),'statevec')).';

end

function grumble(Mhilb)
if ~isnumeric(Mhilb)
    error('Mhilb must be numeric.');
end
end

