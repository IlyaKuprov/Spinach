% Kronecker product between two RCV sparse matrices. Syntax:
%
%                       result=kron(A,B)
%
% Parameters:
%
%    A      - left RCV sparse matrix
%
%    B      - right RCV sparse matrix
%
% Outputs:
%
%    result - RCV sparse matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/kron.m>

function result=kron(A,B)

% Check consistency
grumble(A,B);

% Compute the output dimensions
newRows=A.numRows*B.numRows;
newCols=A.numCols*B.numCols;

% Match location
if (~A.isGPU)&&B.isGPU
    A=gpuArray(A);
end
if (~B.isGPU)&&A.isGPU
    B=gpuArray(B);
end

% Build the Cartesian product of indices and values
[ia,ib]=ndgrid(1:length(A.val),1:length(B.val));
ia=ia(:); ib=ib(:);
newRow=(A.row(ia)-1)*B.numRows+B.row(ib);
newCol=(A.col(ia)-1)*B.numCols+B.col(ib);
newVal=(A.val(ia)).*(B.val(ib));

% Assemble the output object
result=rcv(newCol,newRow,newVal);
result.numRows=newRows;
result.numCols=newCols;
result.isGPU=A.isGPU||B.isGPU;

end

% Consistency enforcement
function grumble(A,B)
if ~isa(A,'rcv')||~isa(B,'rcv')
    error('both inputs must be RCV sparse matrices.');
end
end

% Along the Yangzi River, apes moan ceaselessly.
% My boat has passed ten thousand mounts briskly.
%
% Li Bai

