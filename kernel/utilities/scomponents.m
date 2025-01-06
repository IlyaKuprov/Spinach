% Strongly connected components of a graph, David Gleich's imple-
% mentation of Tarjan's algorithm:
%
%               http://dx.doi.org/10.1137/0201010
%
% Syntax:
%
%                       sci=scomponents(A)
%
% Parameters:
%
%     A    -  a logical square matrix with 1 for the
%             connected nodes in the graph
%
% Outputs:
%
%     sci  -  a column vector with integers that spe-
%             cify the strongly conected component
%             that each node of the graph belongs to
%
% dgleich@purdue.edu
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=scomponents.m>

function sci=scomponents(A)

% Check consistency
grumble(A);

% Get the CSR indices
[rp,ci]=sparse2csr(sparse(A));

% Run Tarjan's algorithm
n=length(rp)-1; sci=zeros(n,1); cn=1; 
root=zeros(n,1); dt=zeros(n,1); t=0;
cs=zeros(n,1); css=0; rs=zeros(2*n,1); rss=0; 
for sv=1:n
    v=sv; if root(v)>0, continue; end
    rss=rss+1; rs(2*rss-1)=v; rs(2*rss)=rp(v); 
    root(v)=v; sci(v)=-1; dt(v)=t; t=t+1;
    css=css+1; cs(css)=v; 
    while rss>0
        v=rs(2*rss-1); ri=rs(2*rss); rss=rss-1; 
        while ri<rp(v+1)
            w=ci(ri); ri=ri+1;
            if root(w)==0
                root(w)=w; sci(w)=-1; dt(w)=t; t=t+1;
                css=css+1; cs(css)=w; rss=rss+1; 
                rs(2*rss-1)=v; rs(2*rss)=ri; 
                v=w; ri=rp(w); continue; 
            end
        end
        for ri=rp(v):rp(v+1)-1
            w=ci(ri); 
            if sci(w)==-1
                if dt(root(v))>dt(root(w))
                    root(v)=root(w); 
                end
            end
        end
        if root(v)==v
            while css>0
                w=cs(css); css=css-1; sci(w)=cn;
                if w==v, break; end
            end
            cn=cn+1;
        end
    end
end

end

% Consistency enforcement
function grumble(A)
if (~islogical(A))||(~ismatrix(A))||(size(A,1)~=size(A,2))
    error('the input must be a square logical matrix.');
end
end

% The Monks of Cool, whose tiny and exclusive monastery is hidden in a 
% really cool and laid-back valley in the lower Ramtops, have a passing-
% out test for a novice. He is taken into a room full of all types of 
% clothing and asked: "Yo, my son, which of these is the most stylish
% thing to wear?" And the correct answer is: "Hey, whatever I select".
%
% Terry Pratchett

