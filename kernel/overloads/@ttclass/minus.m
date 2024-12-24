% Tensor train subtraction operation. Does not perform the actual subtraction
% but instead concatenates the operands until such time as recompression beco-
% mes absolutely necessary. Syntax:
%
%                               c=minus(a,b)
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/minus.m>

function a=minus(a,b)

% Validate the input
if (~isa(a,'ttclass'))||(~isa(b,'ttclass'))
    error('both operands must be tensor train objects.');
elseif (size(a.cores,1)~=size(b.cores,1))||(~all(all(sizes(a)==sizes(b))))
    error('dimension mismatch in tensor trains.');
end

% Write the difference object
a.coeff=[a.coeff -b.coeff];
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

% Meraki (Greek, n.) - the soul, creativity or love put into something;
% the essence of yourself that is put into your work.

