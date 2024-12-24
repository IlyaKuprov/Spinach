% Tensor train addition operation. Does not perform the actual addition, 
% but instead concatenates the operands until such time as recompression 
% becomes absolutely necessary. Syntax:
%
%                              c=plus(a,b)
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/plus.m>

function a=plus(a,b)

% Validate the input
if (~isa(a,'ttclass'))||(~isa(b,'ttclass'))
    error('both operands must be tensor train objects.');
elseif (size(a.cores,1)~=size(b.cores,1))||(~all(all(sizes(a)==sizes(b))))
    error('dimension mismatch in tensor trains.');
end

% Write the sum object
a.coeff=[a.coeff b.coeff];
a.cores=[a.cores b.cores];
a.tolerance=[a.tolerance b.tolerance];

% Filter out zero coeff
pos=find(a.coeff);
if ~isempty(pos)
    a.coeff=a.coeff(pos);
    a.cores=a.cores(:,pos);
    a.tolerance=a.tolerance(pos);
else
    a=0*unit_like(a);
end

end

% Twinkle, twinkle, little star.
% I don't wonder what you are.
% We've got Science, we've got Math!
% You're a giant ball of gas.

