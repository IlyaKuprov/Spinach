% Returns TRUE for non-empty tensor train objects. Syntax:
%
%                answer=isnumeric(tt)
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
% <https://spindynamics.org/wiki/index.php?title=ttclass/isnumeric.m>

function answer=isnumeric(tt)

% Non-empty tensor trains should return true()
if isa(tt,'ttclass')&&(~isempty(tt.cores))
    answer=true();
else
    answer=false();
end

end

% People who think honestly and deeply have a hostile
% attitude towards the public.
%
% Johann Wolfgang von Goethe

