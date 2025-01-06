% Fourier transform axis ticks generator that accounts for the
% periodicity and correctly folds the edge frequency. Syntax:
%
%               ax=ft_axis(offset,sweep,npoints)
%
% Parameters:
%
%    offset   - centre frequency
%
%    sweep    - frequency range
%
%    npoints  - number of points
%
% Outputs:
%
%    ax       - row vector of axis ticks
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ft_axis.m>

function ax=ft_axis(offset,sweep,npoints)

% Check consistency
grumble(offset,sweep,npoints);

% Axis with an extra point
ax=linspace(-sweep/2,sweep/2,npoints+1);

% Odd and even point counts
if mod(npoints,2)==1
    
    % If odd, drop and shift
    ax=ax(2:end)-(ax(2)-ax(1))/2+offset;
    
else
    
    % If even, just drop
    ax=ax(1:(end-1))+offset;
    
end

end

% Consistency enforcement
function grumble(offset,sweep,npoints)
if (~isnumeric(offset))||(~isreal(offset))||...
   (~isscalar(offset))
    error('offset must be a real scalar.');
end
if (~isnumeric(sweep))||(~isreal(sweep))||...
   (~isscalar(sweep))||(sweep<=0)
    error('sweep must be a positive real scalar.');
end
if (~isnumeric(npoints))||(~isreal(npoints))||...
   (~isscalar(npoints))||(mod(npoints,1)~=0)||(npoints<=2)
    error('npoints must be a real integer greater than 2');
end
end

% The role of a conservative thinker is to reassure the people
% that their prejudices are true.
%
% Roger Scruton

