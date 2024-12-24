% Triangular spherical quadrature grids, as per Appendix A.6 of 
% (http://dx.doi.org/10.1016/j.jmr.2014.05.009). Syntax:
%
%        [alps,bets,gams,whts,vorn]=grid_trian(type,n)
%
% Parameters:
%
%      type - 'asg', 'sophe', or 'stoll'
%
%         n - point count parameter
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
% If no outputs are requested, a schematic is drawn.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=grid_trian.m>

function [alps,bets,gams,whts,vorn]=grid_trian(type,n)

% Check consistency
grumble(type,n);

switch type
    
    case 'asg'
        
        % Build the first octant (needs better indexing)
        bets=[]; gams=[];
        for k=0:n
            for j=0:(n-k)
                r=sqrt(k^2+j^2+(n-k-j)^2);
                bets=[bets; acos((n-k-j)/r)]; %#ok<AGROW>
                gams=[gams; atan2(j,k)];      %#ok<AGROW>
            end
        end
        
        % Get point coordinates
        x=sin(bets).*cos(gams);
        y=sin(bets).*sin(gams);
        z=cos(bets);
        
        % Clean up numerics
        mask=isnan(x)|isnan(y)|isnan(z);
        x(mask)=[]; x(abs(x)<sqrt(eps))=0;
        y(mask)=[]; y(abs(y)<sqrt(eps))=0;
        z(mask)=[]; z(abs(z)<sqrt(eps))=0;
        
        % Extend to first quadrant
        z=[z; z(x>0)]; y=[y; y(x>0)]; x=[x; -x(x>0)];
        
        % Extend to top hemisphere
        z=[z; z(y>0)]; x=[x; x(y>0)]; y=[y; -y(y>0)];
        
        % Extend to full sphere
        x=[x; x(z>0)]; y=[y; y(z>0)]; z=[z; -z(z>0)];
        
        % Convert into spherical coordinates
        [gams,elev]=cart2sph(x,y,z); bets=pi/2-elev;
        
        % Alphas are all zero
        alps=zeros(size(gams));
        
    case 'sophe'
        
        % Build the first octant (needs better indexing)
        bets=[]; gams=[];
        for k=0:n
            for j=0:k
                bets=[bets; (pi/2)*(k/n)]; %#ok<AGROW>
                gams=[gams; (pi/2)*(j/k)]; %#ok<AGROW>
            end
        end
        
        % Get point coordinates
        x=sin(bets).*cos(gams);
        y=sin(bets).*sin(gams);
        z=cos(bets);
        
        % Clean up numerics
        mask=isnan(x)|isnan(y)|isnan(z);
        x(mask)=[]; x(abs(x)<sqrt(eps))=0;
        y(mask)=[]; y(abs(y)<sqrt(eps))=0;
        z(mask)=[]; z(abs(z)<sqrt(eps))=0;
        
        % Extend to first quadrant
        z=[z; z(x>0)]; y=[y; y(x>0)]; x=[x; -x(x>0)];
        
        % Extend to top hemisphere
        z=[z; z(y>0)]; x=[x; x(y>0)]; y=[y; -y(y>0)];
        
        % Extend to full sphere
        x=[x; x(z>0)]; y=[y; y(z>0)]; z=[z; -z(z>0)];
        
        % Convert into spherical coordinates
        [gams,elev]=cart2sph(x,y,z); bets=pi/2-elev;
        
        % Cap north and south poles
        bets=[bets; 0; pi]; gams=[gams; 0; 0 ];
        
        % Alphas are all zero
        alps=zeros(size(gams));
        
    case 'stoll'
        
        % Build three SOPHE octants (needs better indexing)
        bets_a=[]; bets_b=[]; bets_c=[];
        gams_a=[]; gams_b=[]; gams_c=[];
        for k=0:n
            for j=0:k % Typo in Stoll's Eq 4.25
                bets_a=[bets_a; (pi/2)*(k/n)];           %#ok<AGROW>
                gams_a=[gams_a; (pi/2)*(j/k)];           %#ok<AGROW>
                bets_b=[bets_b; (pi/2)*((n-j)/n)];       %#ok<AGROW>
                gams_b=[gams_b; (pi/2)*((k-j)/(n-j))];   %#ok<AGROW>
                bets_c=[bets_c; (pi/2)*((n-k+j)/n)];     %#ok<AGROW>
                gams_c=[gams_c; (pi/2)*((n-k)/(n-k+j))]; %#ok<AGROW> 
            end
        end
        
        % Averaging process from Stefan Stoll's thesis
        xa=sin(bets_a).*cos(gams_a); ya=sin(bets_a).*sin(gams_a); za=cos(bets_a);
        xb=sin(bets_b).*cos(gams_b); yb=sin(bets_b).*sin(gams_b); zb=cos(bets_b);
        xc=sin(bets_c).*cos(gams_c); yc=sin(bets_c).*sin(gams_c); zc=cos(bets_c);
        x=xa+yb+zc; y=ya+zb+xc; z=za+xb+yc;
        
        % Clean up numerics
        mask=isnan(x)|isnan(y)|isnan(z);
        x(mask)=[]; x(abs(x)<sqrt(eps))=0;
        y(mask)=[]; y(abs(y)<sqrt(eps))=0;
        z(mask)=[]; z(abs(z)<sqrt(eps))=0;
        
        % Extend to first quadrant
        z=[z; z(x>0)]; y=[y; y(x>0)]; x=[x; -x(x>0)];
        
        % Extend to top hemisphere
        z=[z; z(y>0)]; x=[x; x(y>0)]; y=[y; -y(y>0)];
        
        % Extend to full sphere
        x=[x; x(z>0)]; y=[y; y(z>0)]; z=[z; -z(z>0)];
        
        % Convert into spherical coordinates
        [gams,elev]=cart2sph(x,y,z); bets=pi/2-elev;
        
        % Cap all six poles
        bets=[bets; 0; pi; pi/2; pi/2; pi/2;   pi/2];
        gams=[gams; 0; 0;  0;    pi;   pi/2; 3*pi/2];
        
        % Alphas are all zero
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
function grumble(type,n)
if ~ischar(type)
    error('type must be a character string');
end
if (~isnumeric(n))||(~isreal(n))||...
   (~isscalar(n))||(n<1)||mod(n,1)
    error('n must be a positive real integer');
end
end

% A legend has it that Federico Fellini had once made a bet with Tonio
% Guerra that he would come up with a full-fledged film scenario that 
% would last exactly ten seconds. The following day, the fabled direc-
% tor produced the following sketch: "A woman is watching television.
% A live stream is on from a rocket launch pad, with the countdown: 9,
% 8, 7... Her face shows a storm of emotion. As the last seconds app-
% roach, she picks up a phone and dials a number. As the rocket is le-
% aving the launch pad, she says into the phone: 'Come, he is away!'"

