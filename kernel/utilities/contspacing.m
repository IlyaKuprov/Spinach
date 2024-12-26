% Non-linear adaptive contour spacing. Useful for NMR data where small
% cross-peaks must be adequately contoured next to large diagonal peaks.
% Syntax:
%
%      [all_conts,pos_conts,neg_conts]=...
%                       contspacing(smax,smin,delta,k,signs,ncont)
%
% Parameters:
%
%     smax      - global maximum intensity in the spectrum
%
%     smin      - global minimum intensity in the spectrum
%
%     delta     - minimum and maximum elevation (as a fraction of the
%                 total intensity) of the contours above the baseline.
%                 A good starting value is [0.02 0.2 0.02 0.2]. The
%                 first pair of numbers refers to the positive conto-
%                 urs and the second pair to the negative ones.
%
%     k    - a coefficient that controls the curvature of the contour
%            spacing function: k=1 corresponds to linear spacing and
%            k>1 bends the spacing curve to increase the sampling den-
%            sity near the baseline. A reasonable value is 2.
%
%     signs   - can be set to 'positive', 'negative' or 'both' - this
%               will cause the corresponding contours to be returned.
%
%     ncont   - the number of contours, a reasonable value is 20
%
% Outputs:
%
%     all_conts - all contour levels, a row vector
%
%     pos_conts - positive contour levels, a row vector
%
%     neg_conts - negative contour levels, a row vector
%
% Note: the following functions are used to get contour levels
%
%  pos_conts=delta(2)*smax*linspace(0,1,ncont).^k+smax*delta(1);
%  neg_conts=delta(2)*smin*linspace(0,1,ncont).^k+smin*delta(1);
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=contspacing.m>

function [all_conts,pos_conts,neg_conts]=...
                                 contspacing(smax,smin,delta,k,signs,ncont)

% Check consistency
grumble(smax,smin,delta,k,signs,ncont);

% Compute positive contour levels
if (smax>0)&&(strcmp(signs,'positive')||strcmp(signs,'both'))
    pos_conts=(delta(2)-delta(1))*smax*linspace(0,1,ncont).^k+smax*delta(1);
else
    pos_conts=[]; % Empty if absent or not requested
end

% Compute negative contour levels
if (smin<0)&&(strcmp(signs,'negative')||strcmp(signs,'both'))
    neg_conts=(delta(4)-delta(3))*smin*linspace(0,1,ncont).^k+smin*delta(3);
else
    neg_conts=[]; % Empty if absent or not requested
end

% Merge contour level arrays
all_conts=[neg_conts(end:-1:1) pos_conts];

end

% Consistency enforcement
function grumble(smax,smin,delta,k,signs,ncont)
if (~isnumeric(smax))||(~isscalar(smax))||(~isreal(smax))
    error('smax must be a real scalar.');
end
if (~isnumeric(smin))||(~isscalar(smin))||(~isreal(smin))
    error('smin must be a real scalar.');
end
if (~isnumeric(ncont))||(~isscalar(ncont))||(~isreal(ncont))||...
   (ncont<1)||(mod(ncont,1)~=0)
    error('ncont must be a positive integer.');
end
if (~isnumeric(delta))||(numel(delta)~=4)||(~isreal(delta))||...
    any(delta>1)||any(delta<0)
    error('delta must be a vector with four elements between 0 and 1.');
end
if (~isnumeric(k))||(~isscalar(k))||(~isreal(k))||...
   (k<1)||(mod(k,1)~=0)
    error('k must be a positive integer.');
end
if ~ischar(signs)
    error('signs parameter must be a character string.');
end
end

% "Permit me to point out that you have made 
%  three mistakes in spelling."
%
% Marquis de Favras, after reading his 
% death warrant in 1790

