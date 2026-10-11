% Hash-based row membership for large sparse basis descriptors. Syntax:
%
%                [present,indices]=hash_ismember(A,B)
%
% Parameters:
%
%     A, B     - real, finite double matrices with the same number
%                of columns; either matrix may be sparse or full
%
% Outputs:
%
%     present  - size(A,1)-by-1 logical vector indicating membership
%
%     indices  - size(A,1)-by-1 double vector containing the first
%                matching row in B, or zero when no match exists
%
% Like ismember(A,B,'rows'), rows are compared by numeric value, not
% storage representation. Hashing only nonzero column indices and values
% makes sparse and full rows, including signed zeros, equivalent without
% expanding the basis.
%
% ilya.kuprov@weizmann.ac.il

function [present,indices]=hash_ismember(A,B)

% Check consistency
grumble(A,B);

% Hash the nonzero column indices and values of the query rows
hash_a=repmat(' ',size(A,1),32);
parfor k=1:size(A,1)
    hash_a(k,:)=md5_hash(full(A(k,:)));
end

% Hash the nonzero column indices and values of the reference rows
hash_b=repmat(' ',size(B,1),32);
parfor k=1:size(B,1)
    hash_b(k,:)=md5_hash(full(B(k,:)));
end

% Use hashes to compute membership
[present,indices]=ismember(hash_a,hash_b,'rows');

end

% Consistency enforcement
function grumble(A,B)
if (~isa(A,'double'))||(~ismatrix(A))||(~isreal(A))||any(~isfinite(nonzeros(A)))
    error('A must be a real, finite double matrix.');
end
if (~isa(B,'double'))||(~ismatrix(B))||(~isreal(B))||any(~isfinite(nonzeros(B)))
    error('B must be a real, finite double matrix.');
end
if size(A,2)~=size(B,2)
    error('A and B must have the same number of columns.');
end
end


