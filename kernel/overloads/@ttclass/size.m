% Returns the size of the matrix represented by a tensor train.
% Syntax:
%                          [m,n]=size(tt,dim)
%
% For large spin systems M and N may be too large to fit into 
% the maximum integer permitted by Matlab.
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/size.m>

function varargout=size(tt,dim)

% Multiply up physical dimensions of all cores
if nargin==1
    m=prod(cellfun(@(x)size(x,2),tt.cores),1);
    n=prod(cellfun(@(x)size(x,3),tt.cores),1);
    varargout={[m(1) n(1)]};
elseif dim==1
    m=prod(cellfun(@(x)size(x,2),tt.cores),1);
    varargout={m(1)};
elseif dim==2
    n=prod(cellfun(@(x)size(x,3),tt.cores),1);
    varargout={n(1)};
else
    error('incorrect call syntax.');
end

% Check for infinities
if any(varargout{1}>intmax)
    error('tensor train dimensions exceed Matlab''s intmax.');
end

end

% Will fluorine ever have practical applications? It is very
% difficult to answer this question. I may, however, say in
% all sincerity that I gave this subject little thought when
% I undertook my researches, and I believe that all the chem-
% ists whose attempts preceded mine gave it no more conside-
% ration. A scientific research is a search after truth, and
% it is only after discovery that the question of applicabi-
% lity can be usefully considered.
%
% Henri Moissan

