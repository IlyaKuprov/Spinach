% Generates one descriptor-term operator in XYZ format. Syntax:
%
%                         xyz=term_oper(H,n)
%
% Parameters:
%
%     H   - Hamiltonian action object
%
%     n   - descriptor row number
%
% Outputs:
%
%     xyz - sparse XYZ representation of the descriptor term
%
% ilya.kuprov@weizmann.ac.il
% aditya.dev@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=hamiltonian_action/term_oper.m>

function xyz=term_oper(H,n)

% Check consistency
grumble(H,n);

% Build single-spin operator terms
if H.descr.nS(n)==0

    % Get the single-spin operator
    xyz=operator(H.spin_system,H.descr.opL(n),...
                 {H.descr.nL(n)},H.operator_type,'xyz');

else

    % Get the two-spin operator
    xyz=operator(H.spin_system,[H.descr.opL(n),H.descr.opS(n)],...
                 {H.descr.nL(n),H.descr.nS(n)},H.operator_type,'xyz');

end

end

% Consistency enforcement
function grumble(H,n)
if ~isa(H,'hamiltonian_action')
    error('H must be a Hamiltonian action object.');
end
if (~isnumeric(n))||(~isreal(n))||(~isscalar(n))||...
   (mod(n,1)~=0)||(n<1)||(n>height(H.descr))
    error('n must be a valid descriptor row number.');
end
end

