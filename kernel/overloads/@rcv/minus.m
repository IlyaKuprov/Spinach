% Subtracts one RCV object from another. Syntax:
%
%                    A=minus(A,B)
%
% Parameters:
%
%    A     - left operand
%
%    B     - right operand 
%
% Outputs:
%
%    A     - result A-B as an RCV sparse matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/minus.m>

function A=minus(A,B)

% Check consistency
grumble(A,B);

% Just call plus
A=plus(A,(-1)*B);

end

% Consistency enforcement
function grumble(A,B)
if (~isa(A,'rcv'))&&(~isa(B,'rcv'))
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

