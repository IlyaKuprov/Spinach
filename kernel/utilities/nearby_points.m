% Locates ordered point pairs within an inclusive Euclidean cutoff. Syntax:
%
%             [from,to]=nearby_points(xyz,images,r0)
%
% Parameters:
%
%    xyz, images - arrays of Cartesian coordinates, one point per row
%             r0 - positive cutoff in the same length units
%
% Outputs:
%
%       from, to - column vectors of row indices into xyz and images,
%                  respectively; coincident points are included
%
% ilya.kuprov@weizmann.ac.il
% talos@spindynamics.org

function [from,to]=nearby_points(xyz,images,r0)

% Check consistency
grumble(xyz,images,r0);

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
    offset=offsets(m,:);
    target=origins+offset;
    [present,which_bin]=ismember(target,bins,'rows');
    atoms=find(present&all((offset==0)|(target~=origins),2));
    block=cell(numel(atoms),1);
    for n=1:numel(atoms)

        % Test Euclidean distances against the cutoff
        a=atoms(n);
        b=members{which_bin(a)};
        dv=images(b,:)-xyz(a,:);
        keep=vecnorm(dv,2,2)<=r0;
        block{n}=[repmat(a,nnz(keep),1) b(keep)];

    end
    hits{m}=vertcat(block{:});

end

% Concatenate the occupied pairs, including the empty case
hits=vertcat(zeros(0,2),hits{:});
from=hits(:,1);
to=hits(:,2);

end

% Consistency enforcement
function grumble(xyz,images,r0)
if (~isnumeric(xyz))||(~isreal(xyz))||(size(xyz,2)~=3)||...
   (~isnumeric(images))||(~isreal(images))||(size(images,2)~=3)
    error('coordinates must be real matrices with three columns.');
end
if (~isnumeric(r0))||(~isreal(r0))||(~isscalar(r0))||(r0<=0)
    error('cutoff must be a positive real scalar.');
end
end


