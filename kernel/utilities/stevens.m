% Extended Stevens operators. Syntax:
%
%                      S=stevens(mult,k,q)
%
% Parameters:
%
%   mult - multiplicity of the spin in question
%
%   k    - Stevens operator rank, a non-negative
%          integer
%
%   q    - Stevens operator projection, an
%          integer between -k and k
%
% Outputs:
%
%   S    - Stevens operator matrix
%
% Note: for historical reasons, the definition of Stevens operators
%       is irregular and must rely on explicitly stockpiled coeffi-
%       cients. Only ranks smaller or equal to 12 are available.
%
% ilya.kuprov@weizmann.ac.il
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=stevens.m>

function S=stevens(mult,k,q)

% Check consistency
grumble(mult,k,q);

% Stockpile the ugly integer coefficients (there is no systematic formula)
C{13}=[1916006400 958003200 958003200 31933440 3991680 1995840 31680 15840 1584 264 24 12 1 ];
C{12}=[319334400  79833600  79833600  13305600 2661120 23760   7920  1320  1320 22  22 1 ];
C{11}=[14515200   7257600   1209600   604800   86400   2880    360   180   20   10  1 ];
C{10}=[1451520    725700    725700    60480    60480   864     288   18    18   1 ];
C{9} =[80640      40320     40320     6720     672     336     16    8     1 ];
C{8} =[40320      5040      1680      168      168     14      14    1 ];
C{7} =[2880       1440      360       60       12      6       1 ];
C{6} =[480        240       240       10       10      1 ];
C{5} =[48         24        8         4        1 ];
C{4} =[24         6         6         1 ];
C{3} =[4          2         1 ];
C{2} =[2          1 ];
C{1} =[1 ]; %#ok<*NBRAK>

% Get Pauli matrices
L=pauli(mult);

% Get the top state
S=L.p^k;

% Commute down
for n=abs(q):(k-1)
    S=L.m*S-S*L.m;
end 

% Distinguish odd and even indices
if (~logical(mod(k,2)))&&(logical(mod(q,2)))
    coeff=((-1)^(k-q))/C{k+1}(abs(q)+1)/2;
else
    coeff=((-1)^(k-q))/C{k+1}(abs(q)+1);
end

% Distinguish positive and negative indices
if (q>0)||(q==0)
    S=coeff*(S+S')/2;
else 
    S=coeff*(S-S')/2i;
end

end

% Consistency enforcement
function grumble(mult,k,q)
if (~isnumeric(mult))||(~isreal(mult))||...
   (~isscalar(mult))||(mod(mult,1)~=0)||(mult<0)
    error('the multiplicity must be a positive real integer.');
end
if (~isnumeric(k))||(~isreal(k))||(~isfinite(k))||...
   (~isscalar(k))||(mod(k,1)~=0)||(k<0)||(k>12)
    error('k must be a real integer between 0 and 12.');
end
if (~isnumeric(q))||(~isreal(q))||(~isfinite(q))||...
   (~isscalar(q))||(mod(q,1)~=0)||(q<-k)||(q>k)
    error('q must be an integer between -k and k.');
end
end

% When nine hundred years old you reach, look as 
% good you will not.
%
% Master Yoda

