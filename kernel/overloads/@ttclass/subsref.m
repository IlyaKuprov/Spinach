% Dot and bracket property specifications for the tensor train class.
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.uk
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
            ind=ttclass_ind2sub(siz(:,1),reference(1).subs{1});
            jnd=ttclass_ind2sub(siz(:,2),reference(1).subs{2});
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

% #NGRUM #NHEAD