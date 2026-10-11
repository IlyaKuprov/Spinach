% Hash-based row membership for large sparse basis descriptors. Syntax:
%
%                [present,indices]=hash_ismember(A,hash_b)
%
% Parameters:
%
%     A        - real, finite double query matrix, sparse or full
%
%     hash_b   - N-by-32 character array of reference row MD5 hashes,
%                each computed as md5_hash(full(B(k,:))) for a real,
%                finite double reference matrix B with size(A,2) columns
%
% Outputs:
%
%     present  - size(A,1)-by-1 logical vector indicating membership
%
%     indices  - size(A,1)-by-1 double vector containing the first
%                matching row in B, or zero when no match exists
%
% Full rows are hashed, independently of sparse storage. Reference row
% hashes are supplied by the caller and are not recomputed. MD5 collisions
% are not checked; column counts cannot be checked from hashes alone.
%
% ilya.kuprov@weizmann.ac.il

function [present,indices]=hash_ismember(A,hash_b)

% Check consistency
grumble(A,hash_b);

% Hash full query rows
hash_a=repmat(' ',size(A,1),32);
parfor k=1:size(A,1)
    hash_a(k,:)=md5_hash(full(A(k,:)));
end

% Use hashes to compute membership
[present,indices]=ismember(hash_a,hash_b,'rows');

end

% Consistency enforcement
function grumble(A,hash_b)
if (~isa(A,'double'))||(~ismatrix(A))||(~isreal(A))||any(~isfinite(nonzeros(A)))
    error('A must be a real, finite double matrix.');
end
if (~ischar(hash_b))||(~ismatrix(hash_b))||(size(hash_b,2)~=32)
    error('hash_b must be an N-by-32 character array of reference row hashes.');
end
end


