% Returns TRUE for real-valued tensor train objects. Syntax:
%
%                answer=isreal(tt)
%
% Parameters:
%
%    tt - tensor train object
%
% Outputs:
%
%    answer - logical true when all coefficients and core
%             elements of the tensor train are real
%
% d.savostyanov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/isreal.m>

function answer=isreal(tt)

% Non-empty tensor trains should return true()
if isa(tt,'ttclass')

    % Check coefficient first 
    answer=all(isreal(tt.coeff));
    
    % If the coefficients are real, check the cores
    if answer
       for n=1:tt.ntrains
           for k=1:tt.ncores
               answer=isreal(tt.cores{k,n});
               if ~answer, return; end
           end
       end 
    end

else
    
    % Complain and bomb out
    error('input is not a ttclass.');

end

end

% Democracy is a pathetic belief in the collective wisdom
% of individual ignorance.
%
% H.L. Mencken

% #NGRUM

