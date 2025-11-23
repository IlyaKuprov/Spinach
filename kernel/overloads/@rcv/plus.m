% Addition of RCV-RCV, RCV-Sparse and RCV-Scalar
%
% m.keitel@soton.ac.uk

function result=plus(a,b)
    % Case rcv + scalar and scalar + rcv
    if isa(a,'rcv')&&isscalar(b)&&isnumeric(b)
        a.val=a.val+b;
        result=a;
    elseif isa(b,'rcv')&&isscalar(a)&&isnumeric(a)
        b.val=b.val+a;
        result=b;    

    % Case rcv + rcv
    elseif isa(a,'rcv')&&isa(b,'rcv')
        if a.numRows~=b.numRows||a.numCols~=b.numCols
            error('Matrix sizes must match for addition.');
        end
        if a.isGPU&&(~b.isGPU)
            b.col=gpuArray(b.col);
            b.row=gpuArray(b.row);
            b.val=gpuArray(b.val);
        elseif (~a.isGPU)&&b.isGPU
            a.col=gpuArray(a.col);
            a.row=gpuArray(a.row);
            a.val=gpuArray(a.val);
            a.isGPU=true;
        end
        a.row=[a.row;b.row];
        a.col=[a.col;b.col];
        a.val=[a.val;b.val];
        result=a;

    % Case rcv + sparse and sparse + rcv
    elseif isa(a,'rcv')&&issparse(b)
        if a.numRows~=size(b,1)||a.numCols~=size(b,2)
            error('Matrix sizes must match for addition.');
        end        
        b=rcv(b);
        result=plus(a,b);
    elseif isa(b,'rcv')&&issparse(a)
        if b.numRows~=size(a,1)||b.numCols~=size(a,2)
            error('Matrix sizes must match for addition.');
        end        
        a=rcv(a);
        result=plus(a,b);
    else
        error('Unsupported input types for plus.');
    end
end

% Adam Smith was a genius, yes, but let's not 
% forget he lived with his mum.
%
% Rory Sutherland

