% Returns true for unit matrices. The test is designed to be
% computationally affordable. Syntax:
%
%                       verdict=iseye(M)
%
% Parameters:
%
%    M        - a matrix
%
% Outputs:
%
%    verdict  - true or false
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=iseye.m>

function verdict=iseye(M)

% Check consistency
grumble(M);

% Run the checks
if size(M,1)~=size(M,2)
    
    % Not even square
    verdict=false();
    
elseif ~isdiag(M)
    
    % Not even diagonal
    verdict=false();
    
else
    
    % Test vector
    a=randn(size(M,2),1);
    
    % Compare with unit
    if nnz(M*a-a)~=0
        
        % Test failed
        verdict=false();
        
    else
        
        % Actually unit
        verdict=true();
        
    end
    
end

end

% Consistency enforcement
function grumble(M)
if ~isnumeric(M)
    error('M must be numeric.');
end
end

% Some candidates also reproduced an unnecessary derivation of 
% the quantum Hall effect, the question having clearly set off 
% some sort of Pavlovian response.
%
% Oxford Chemistry Examiners' Report 2015

