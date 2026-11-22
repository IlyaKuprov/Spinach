% Finds ordered pairs of distinct atoms strictly within a distance cutoff.
% Syntax:
%
%                 pairs=conmat(xyz,images,r0)
%
% Parameters:
%
%    xyz, images - N-by-3 coordinates of the atoms in the reference
%                  cell and one translated cell, in the same order
%
%             r0 - non-negative cutoff, in coordinate units
%
% Outputs:
%
%          pairs - cell rows {from,to,displacement}, from the reference
%                  atom to the image atom, ordered by to, then by from
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=conmat.m>

function pairs=conmat(xyz,images,r0)

% Check consistency
grumble(xyz,images,r0);

% No points can be strictly within a zero cutoff
pairs=cell(0,3);
if isempty(xyz)||(r0==0), return; end

% Group image atoms by occupied cutoff-sized cell
[bins,~,group]=unique(floor(images/r0),'rows');
members=accumarray(group,(1:size(images,1))',[],@(x){x});
origins=floor(xyz/r0);

% Search neighbouring cells for each atom
[x,y,z]=ndgrid(-1:1);
offsets=[x(:) y(:) z(:)];
hits=cell(size(offsets,1),1);
for m=1:size(offsets,1)

    % Locate occupied image cells at this offset
    [present,which_bin]=ismember(origins+offsets(m,:),bins,'rows');
    atoms=find(present);
    block=cell(numel(atoms),1);
    for n=1:numel(atoms)

        % Test Euclidean distances and exclude all self images
        a=atoms(n);
        b=members{which_bin(a)};
        dv=images(b,:)-xyz(a,:);
        keep=(a~=b)&(sqrt(sum(dv.^2,2))<r0);
        block{n}=[repmat(a,nnz(keep),1) b(keep) dv(keep,:)];

    end
    hits{m}=vertcat(block{:});

end

% Gather ordered pairs and their displacement vectors
hits=vertcat(hits{:});
if ~isempty(hits)
    hits=sortrows(hits,[2 1]);
    pairs=[num2cell(hits(:,1:2)) num2cell(hits(:,3:5),2)];
end

end

% Consistency enforcement
function grumble(xyz,images,r0)
if (~isnumeric(xyz))||(~isreal(xyz))||(~ismatrix(xyz))||(size(xyz,2)~=3)
    error('xyz must be a real matrix with three columns.');
end
if (~isnumeric(images))||(~isreal(images))||(~isequal(size(images),size(xyz)))
    error('images must be a real matrix of the same size as xyz.');
end
if (~isnumeric(r0))||(~isreal(r0))||(~isscalar(r0))||(~(r0>=0))
    error('r0 must be a non-negative real scalar.');
end
end

% It's always good to be underestimated.
%
% Donald Trump


