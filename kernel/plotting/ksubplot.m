% House style settings for Matlab figures; a product of much
% experience with academic publication aesthetics. This is a
% replacement for Matlab's subplot function with narrower pa-
% nel and figure margins. Syntax:
%
%                      ksubplot(m,n,k)
%
% Parameters:
%
%    m,n,k - same arguments as for Matlab's subplot()
%
% Outputs:
% 
%    creates or updates the subplot layout 
%    of the current figure
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ksubplot.m>

function ksubplot(m,n,k)

% Check consistency
grumble(m,n,k);

% Gap and margin settings
gap_vert=0.01; gap_horz=0.01;
marg_b=0.05; marg_t=0.05;
marg_l=0.05; marg_r=0.05;

% Linear index to row, col
[col,row]=ind2sub([n,m],k);  

% Fractional dimensions of the subplot
height=(1-(marg_b+marg_t)-(m-1)*gap_vert)/m;
width =(1-(marg_l+marg_r)-(n-1)*gap_horz)/n;

% Fractional position of the subplot
bottom=(m-row)*(height+gap_vert)+marg_b;
left  =(col-1)*(width +gap_horz)+marg_l;

% Matlab subplot function call
subplot('Position',[left bottom width height]);

end

% Consistency enforcement
function grumble(m,n,k)
if (~isnumeric(m))||(~isreal(m))||(~isscalar(m))||(mod(m,1)~=0)||(m<1)||...
   (~isnumeric(n))||(~isreal(n))||(~isscalar(n))||(mod(n,1)~=0)||(n<1)||...
   (~isnumeric(k))||(~isreal(k))||(~isscalar(k))||(mod(k,1)~=0)||(k<1)
    error('m, n, k must be positive real integers.');
end
if k>m*n, error('k cannot exceed m*n.'); end
end

% I have never seen a thin person 
% drinking Diet Coke.
%
% Donald Trump

