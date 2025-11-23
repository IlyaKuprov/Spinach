% Subtracts one RCV object from another. Syntax:
%
%                    a=minus(a,b)
%
% Parameters:
%
%    a     - left operand
%
%    b     - right operand 
%
% Outputs:
%
%    a     - result of a-b
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/minus.m>

function a=minus(a,b)

% Check consistency
grumble(a,b);

% Just call plus
a=plus(a,(-1)*b);

end

% Consistency enforcement
function grumble(a,b)
if (~isa(a,'rcv'))&&(~isa(b,'rcv'))
    error('at least one input must be an RCV sparse matrix.');
end
end

% "My cat had been suffering from severe illness over the past month
%  or so. This had meant that he had needed increasingly hands-on
%  care. Due to a terminal diagnosis the decision to put him to
%  sleep was made; that took place on Monday 11th April.
%
%  This has caused me grievance, affected my moods and my general
%  ability to concentrate and perform my academic tasks. I am un-
%  able to provide suporting evidence as I do not have a receipt
%  for the final vet appointment, and, as of yet, have no death
%  certificate or other evidence of death. I can provide photos
%  of the cat."
%
% A 2022 special considerations request by a
% Southampton University student, requesting
% more lenient assessment.

