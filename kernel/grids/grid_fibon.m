% Fibonacci type spherical quadrature grids, as per Appendix 
% A.5 of http://dx.doi.org/10.1016/j.jmr.2014.05.009 Syntax:
%
%      [alps,bets,gams,whts,vorn]=grid_fibon(type,parm)
%
% Parameters:
%
%      type - 'fibonacci', 'zcw', or 'zcwn'
%
%      parm - point count parameter, the resulting
%             grid will have 2*n+1 points ('fib'),
%             fibonacci(n+2) points ('zcw'), or n
%             points (zcwn).
%
% Outputs:
%
%      alps - alpha Euler angles of the grid (radians),
%             zeros because these are two-angle grids
%
%      bets - beta Euler angles of the grid (radians)
%
%      gams - gamma Euler angles of the grid (radians)
% 
%      whts - Voronoi tessellation body angle weights
%
%      vorn - a cell array of matrices containing the
%             coordinates of the vertices of the Voro-
%             noi polyhedra
%
% Note: if no outputs are requested, a schematic is drawn.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=grid_fibon.m> 

function [alps,bets,gams,whts,vorn]=grid_fibon(type,parm)

% Check consistency
grumble(type,parm);

switch type
    
    case 'fib'
        
        % Golden ratio
        phi=(1+sqrt(5))/2;
        
        % Point index
        k=(-parm:parm)';
        
        % Fibonacci angles
        bets=acos(2*k/(2*parm+1));
        gams=2*pi*k/phi;
        alps=zeros(size(gams));
        
    case 'zcw'
        
        % Fibonacci numbers
        Fmp2=fibonacci(parm+2);
        Fm=fibonacci(parm); 
        
        % Point index
        k=(0:(Fmp2-1))';
        
        % Zaremba's grid
        bets=acos(2*k/Fmp2-1);
        gams=2*pi*k*Fm/Fmp2;
        alps=zeros(size(gams));
    
    case 'zcwn'
        
        % Golden ratio
        phi=(1+sqrt(5))/2;

        % Craciun's grid
        bets=acos(2*((0:(parm-1))'/parm)-1);
        gams=2*pi*(0:(parm-1))'/phi^2;
        alps=zeros(size(gams));
        
    otherwise
        
        % Complain and bomb out
        error('unknown grid type');
        
end

% Weights are expensive
if (nargout>3)||(nargout==0)
    
    % Get point coordinates
    x=sin(bets).*cos(gams);
    y=sin(bets).*sin(gams);
    z=cos(bets);
    
    % Run Voronoi tessellation
    [~,~,vorn,whts]=voronoisphere([x'; y'; z']);
    
    % Get weights
    whts=whts/(4*pi);
    
end

% If no output requested, plot the grid
if nargout==0, grid_plot(x,y,z,vorn); end

end

% Consistency enforcement
function grumble(type,parm)
if ~ischar(type)
    error('type must be a character string');
end
if (~isnumeric(parm))||(~isreal(parm))||...
   (~isscalar(parm))||(parm<1)||mod(parm,1)
    error('parm must be a positive real integer');
end
end

% On the whole, human beings want to be good, but not too 
% good and not quite all of the time.
%
% George Orwell

