% Yamaguchi equation estimate of exchange coupling from
% a broken-symmetry DFT calculation on a bistrityl bira-
% dical with an alkinyl linker from Olav Schiemann.
%
% i.kuprov@soton.ac.uk

function yamaguchi()

% Read Gaussian logs
props_sing=gparse('biradical_singlet.log');
props_trip=gparse('biradical_triplet.log');

% Call Yamaguchi equation
J=brokensymm(props_sing,props_trip);
disp(['Exchange coupling: ' num2str(J/1e9) ' GHz']);

end

