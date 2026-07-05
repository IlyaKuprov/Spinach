% Sparse matrix row-sorting permutation utility. Syntax:
%
%                    idx=spsortrows(A)
%
% Parameters:
%
%      A   - sparse real double matrix
%
% Outputs:
%
%      idx - row permutation index, matching the second
%            output of Matlab's sortrows(A)
%
% This file is a Matlab fallback for the compiled MEX function.
%
% ilya.kuprov@weizmann.ac.il

function idx=spsortrows(A)

% Check consistency
grumble(A);

% Return Matlab reference permutation
[~,idx]=sortrows(A);

end

% Consistency enforcement
function grumble(A)
if (~isnumeric(A))||(~issparse(A))||(~isreal(A))||...
   (~isa(A,'double'))||(~ismatrix(A))
    error('A must be a sparse real double matrix.');
end
end

% Nihilistic Password Security Questions
%
%   What is the name of your least favorite child?
%
%   In what year did you abandon your dreams?
%
%   What is the maiden name of your father's mistress?
% 
%   At what age did your childhood pet run away?
%
%   What was the name of your favorite unpaid internship?
%
%   In what city did you first experience ennui?
%
%   What is your ex-wife's newest last name?
%
%   What sports team do you fetishise to avoid meaningful discussion with others?
%
%   What is the name of your favorite cancelled TV show?
%
%   What was the middle name of your first rebound?
%
%   On what street did you lose your childlike sense of wonder?
%
%   When did you stop trying?
%
% Soheil Rezayazdi

