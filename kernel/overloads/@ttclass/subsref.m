% Dot and bracket property specifications for the tensor train class.
% Syntax:
%
%                answer=subsref(ttrain,reference)
%
% Parameters:
%
%    ttrain    - tensor train object
%
%    reference - Matlab subscript reference structure
%
% Outputs:
%
%    answer    - requested tensor train property, scalar
%                matrix element, or nested subscript result
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/subsref.m>

function answer=subsref(ttrain,reference)

switch reference(1).type
    
    % Methods and properties
    case '.'
        
        % Return the output requested
        switch reference(1).subs 
            case 'ncores',    answer=ttrain.ncores;
            case 'ntrains',   answer=ttrain.ntrains;
            case 'sizes',     answer=ttrain.sizes; 
            case 'ranks',     answer=ttrain.ranks;           
            case 'coeff',     answer=ttrain.coeff;
            case 'cores',     answer=ttrain.cores;
            case 'tolerance', answer=ttrain.tolerance;
            otherwise,        error(['unknown field reference ',reference(1).subs]);
        end
    
    % Matrix element extraction
    case '()'
        
        % Start with zero
        answer=0;
        
        % Convert indices
        if numel(reference(1).subs)~=2
            error('exactly two indices required to evaluate tt(j,k) element');
        elseif isscalar(reference(1).subs{1})&&isscalar(reference(1).subs{2})
            siz=sizes(ttrain);
            row_ind=reference(1).subs{1}; col_ind=reference(1).subs{2};
            if islogical(row_ind), row_ind=double(row_ind); end
            if islogical(col_ind), col_ind=double(col_ind); end
            grumble(row_ind,col_ind,siz);
            ind=ttclass_ind2sub(siz(:,1),row_ind);
            jnd=ttclass_ind2sub(siz(:,2),col_ind);
        else
            error('advanced indexing is not implemented for tensor trains.');
        end
        
        % Multiply up the tensor train
        [ncores,ntrains]=size(ttrain.cores);
        for n=1:ntrains
            x=ttrain.cores{ncores,n}(:,ind(ncores),jnd(ncores),:);
            for k=(ncores-1):(-1):1
                x=ttrain.cores{k,n}(:,ind(k),jnd(k),:)*x;
            end
            answer=answer+ttrain.coeff(n)*x(1,1);
        end
        
    otherwise
        
        % Complain and bomb out
        error('unknown subscript reference type.');
        
end

% Allow nested indexing
if numel(reference)>1
    answer=subsref(answer,reference(2:end));
end

end

% Consistency enforcement
function grumble(row_ind,col_ind,siz)
if (~isnumeric(row_ind))||(~isreal(row_ind))||...
   (mod(row_ind,1)~=0)||(row_ind<1)||(row_ind>flintmax)
    error('row index must be a positive integer not exceeding flintmax.');
end
if (~isnumeric(col_ind))||(~isreal(col_ind))||...
   (mod(col_ind,1)~=0)||(col_ind<1)||(col_ind>flintmax)
    error('column index must be a positive integer not exceeding flintmax.');
end
leftover=row_ind-1;
for n=1:size(siz,1)
    leftover=floor(leftover/siz(n,1));
end
if leftover>0
    error('row index must be a positive integer within matrix dimensions.');
end
leftover=col_ind-1;
for n=1:size(siz,1)
    leftover=floor(leftover/siz(n,2));
end
if leftover>0
    error('column index must be a positive integer within matrix dimensions.');
end
end

% Index to subscript transformation
function ivec=ttclass_ind2sub(siz,ind)
d=numel(siz); ind=ind-1; ivec=zeros(d,1);
for k=d:-1:1
    ivec(k)=mod(ind,siz(k));
    ind=floor(ind/siz(k));
end
ivec=ivec+1;
end

% "A cactus is a very disappointed cucumber."
%
% A Russian saying


