% Forward linear prediction. Syntax:
%
%           y=lpredict(x,npcoeffs,npredps)
%
% Parameters:
%
%      x        - input data, a column vector
%
%      npcoeffs - number of predictor coefficients,
%                 must be greater than 1
%
%      npredps  - number of data points to predict
%
% Outputs:
%
%      y        - predicted points
%
% malcolm.lidierth@kcl.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=lpredict.m>

function y=lpredict(x,npcoeffs,npredps)

% Check consistency
grumble(x,npcoeffs,npredps);

% Store and subtract the mean
mean_x=mean(x); x=x-mean_x;

% Store and scale by stdev
std_x=std(x); x=x/std_x;

% Get linear predictor coefficients
a=lpc(x,npcoeffs); lpcoeffs=-a(2:end);

% Pre-allocate output
y=zeros(npredps,1);

% First value only uses x
y(1)=lpcoeffs*x(end:-1:(end-npcoeffs+1));

% Then some x and some y
for k=2:min(npcoeffs,npredps)
    y(k)=lpcoeffs*[y((k-1):-1:1); ...
                   x(end:-1:(end-npcoeffs+k))];
end

% Then keep going using y
for k=(npcoeffs+1):npredps
    y(k)=lpcoeffs*y(k-1:-1:k-npcoeffs);
end

% Scale and shift back
y=std_x*y; y=y+mean_x;

end

% Consistency enforcement
function grumble(x,npcoeffs,npredps)
if (~isnumeric(x))||(~isreal(x))||(~iscolumn(x))
    error('x must be a real column vector.');
end
if (~isnumeric(npcoeffs))||(~isreal(npcoeffs))||...
   (~isscalar(npcoeffs))||(npcoeffs<2)||(mod(npcoeffs,1)~=0)
    error('npcoeffs must be a real integer greater than 1.');
end
if (~isnumeric(npredps))||(~isreal(npredps))||...
   (~isscalar(npredps))||(npredps<1)||(mod(npredps,1)~=0)
    error('npredps must be a positive real integer.');
end
end

% There is a serious tendency towards capitalism
% among the well-to-do peasants.
%
% Mao Zedong

