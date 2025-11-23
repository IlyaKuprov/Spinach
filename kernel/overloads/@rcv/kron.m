% Kronecker Product of RCV objects
%
% m.keitel@soton.ac.uk

function result=kron(A,B)

    if ~isa(A,'rcv')||~isa(B,'rcv')
        error('Both inputs must be CRV objects');
    end

    % Compute new dimensions
    newRows=A.numRows*B.numRows;
    newCols=A.numCols*B.numCols;

    % Cartesian product of row/col indices and value multiplication
    % Ensure everything is on GPU
    if ~A.isGPU&&B.isGPU
        A=gpuArray(A);
    end
    if ~B.isGPU&&A.isGPU
        B=gpuArray(B);
    end
    
    [ia,ib]=ndgrid(1:length(A.val),1:length(B.val));
    ia=ia(:); ib=ib(:);

    newRow=(A.row(ia)-1)*B.numRows+B.row(ib);
    newCol=(A.col(ia)-1)*B.numCols+B.col(ib);
    newVal=(A.val(ia)).*(B.val(ib));

    result=rcv(newCol,newRow,newVal);
    result.numRows=newRows;
    result.numCols=newCols;
    result.isGPU=A.isGPU||B.isGPU;
    
end

% Along the Yangzi River, apes moan ceaselessly.
% My boat has passed ten thousand mounts briskly.
%
% Li Bai

