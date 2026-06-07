% Returns TRUE for non-empty tensor train objects. Syntax:
%
%                answer=ismatrix(tt)
%
% Parameters:
%
%    tt - tensor train object
%
% Outputs:
%
%    answer - logical true for non-empty tensor train objects
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/ismatrix.m>

function answer=ismatrix(tt)

% Non-empty tensor trains should return true()
if isa(tt,'ttclass')&&(~isempty(tt.cores))
    answer=true();
else
    answer=false();
end

end

% The stronger the house, the greater the immigration.
%
% The Law of Three Little Pigs

% #NGRUM

