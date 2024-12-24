% Generates a tensor train representation of a matrix with random
% tensor train cores, same physical index topology at the tensor
% train supplied, and user-specified bond ranks. Syntax:
%
%                        tt=rand(tt,ttrank)
%
% Parameters:
%
%    tt     - a tensor train object
%
%    ttrank - bond rank, a positive integer
%
% Outputs:
%
%    tt     - a tensor train object
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/rand.m>

function tt=rand(tt,ttrank)

% Check consistency
grumble(tt,ttrank);

% Read tensor train sizes
[d,~]=size(tt.cores); sz=sizes(tt);

% Reallocate cores
tt.cores=cell(d,1);

% Fill the cores with random elements
tt.cores{1,1}=rand(1,sz(1,1),sz(1,2),ttrank);
for k=2:d-1
    tt.cores{k,1}=rand(ttrank,sz(k,1),sz(k,2),ttrank);
end
tt.cores{d,1}=rand(ttrank,sz(d,1),sz(d,2),1);

% Unit coefficient and zero tolerance
tt.coeff=1; tt.tolerance=0;

end

% Consistency enforcement
function grumble(tt,ttrank)
if ~isa(tt,'ttclass')
    error('tt must be a tensor train object.');
end
if (~isnumeric(ttrank))||(~isreal(ttrank))||...
   (~isscalar(ttrank))||(ttrank<1)||(mod(ttrank,1)~=0)
    error('ttrank must be a positive real integer.');
end
end

% The problem with inviting professors to dinner parties 
% is that they're used to talking for a minimum of 50 mi-
% nutes at a time.
%
% Anonymous philosophy professor

