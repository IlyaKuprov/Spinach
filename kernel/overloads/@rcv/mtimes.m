% Matrix multiplication
%
% m.keitel@soton.ac.uk

function result=mtimes(a,b)
    % Case: scalar * rcv or rcv * scalar
    if isa(a,'rcv')&&isscalar(b)&&isnumeric(b)
        a.val=a.val*b;
        result=a;

    elseif isa(b,'rcv')&&isscalar(a)&&isnumeric(a)
        b.val=b.val*a;
        result=b;

    % Case: rcv * rcv
    elseif isa(a,'rcv')&&isa(b,'rcv')
        if a.numCols~=b.numRows
            error('Inner matrix dimensions must agree.');
        end
        result=sparse(a)*sparse(b);

    % Case: rcv * sparse
    elseif isa(a,'rcv')&&issparse(b)
        if a.numCols~=size(b,1)
            error('Inner matrix dimensions must agree.');
        end
        result=sparse(a)*b;
    elseif issparse(a)&&isa(b,'rcv')
        if size(a,2)~=b.numRows
            error('Inner matrix dimensions must agree.');
        end
        result=a*sparse(b);
    else
        error('Unsupported input types for mtimes.');
    end

end

% When the Earl of Sandwich, speaking in parliament, told John Wilkes that
% the latter would either die on the gallows or of the pox, Wilkes politely
% responded that it would depend on whether he embraced the earl's princi-
% ples or his mistress.
%
% Taki Theodoracopoulos

