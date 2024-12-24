% Calculates finite difference weights for numerical derivatives, 
% including order 0, which amounts to interpolation. Syntax:
%
%         w=fdweights(target_point,grid_points,max_order)
%
% Parameters:
%
%     target_point - the point at which the derivative 
%                    is required
%
%     grid_points  - the points at which the function 
%                    is given
%
%     max_order    - maximum derivative order
%
% Outputs:
%
%     w            - finite difference coefficient array
%                    with the coefficients for the succes-
%                    sive derivatives in rows
%
% fornberg@colorado.edu
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=fdweights.m>

function w=fdweights(target_point,grid_points,max_order)

% Check consistency
grumble(target_point,grid_points,max_order);

% Compute the weights
n=length(grid_points); w=zeros(max_order+1,n); 
w1=1; w4=grid_points(1)-target_point; w(1,1)=1;
for i=2:n
   mn=min(i,max_order+1); w2=1; w5=w4; w4=grid_points(i)-target_point;
   for j=1:(i-1)
      w3=grid_points(i)-grid_points(j); w2=w2*w3;
      if j==(i-1) 
         w(2:mn,i)=w1*((1:(mn-1))'.*w(1:(mn-1),i-1)-w5*w(2:mn,i-1))/w2;
         w(1,i)=-w1*w5*w(1,i-1)/w2;
      end
      w(2:mn,j)=(w4*w(2:mn,j)-(1:(mn-1))'.*w(1:(mn-1),j))/w3;
      w(1,j)=w4*w(1,j)/w3;
   end
   w1=w2;
end

end

% Consistency enforcement
function grumble(target_point,grid_points,max_order)
if (~isnumeric(target_point))||(~isnumeric(grid_points))||(~isnumeric(max_order))||...
   (~isreal(target_point))||(~isreal(grid_points))||(~isreal(max_order))
    error('all input arguments must be numeric and real.');
end
if numel(target_point)~=1
    error('target_point must be a scalar.');
end
if ~isvector(grid_points)
    error('grid_points must be a vector.');
end
if (target_point>max(grid_points))||(target_point<min(grid_points))
    error('the target point should be inside the grid.');
end
if ~issorted(grid_points)
    error('grid_points array must be sorted in ascending order.');
end
if mod(max_order,1)~=0
    error('max_order must be an integer.');
end
if max_order>=numel(grid_points)
    error('max_order must be smaller than the number of grid points.');
end
end

% Conscience and cowardice are really the same things. Conscience
% is the trade name of the firm.
%
% Oscar Wilde

