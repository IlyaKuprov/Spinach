% Converts the Weblab one-cone model parameters (see weblab_cone.png)
% into NQI tensors used by Spinach. Syntax:
%
%     [Q1,Q2]=weblab2nqi(C_q,eta_q,I,alpha,theta,phi)
%     [Q1,Q2,Q3]=weblab2nqi(C_q,eta_q,I,alpha,theta,phi)
%     [Q1,Q2,Q3,Q4]=weblab2nqi(C_q,eta_q,I,alpha,theta)
%     [Q1,Q2,Q3,Q4,Q5,Q6]=weblab2nqi(C_q,eta_q,I,alpha,theta)
%
% Parameters:
%
%     C_q          - quadrupolar coupling constant e^2*q*Q/h
%                    in Hz
%
%     eta_q        - quadrupolar tensor asymmetry parameter
%
%     I            - spin quantum number 
%
%     alpha
%     theta
%     phi          - the three angles of Weblab cone model
%                    (see weblab_cone.png), in radians; the
%                    four- and six-site modes place the sites
%                    on fixed azimuth grids, [0 1 2 3]*pi/2
%                    and [0 1 2 3 4 5]*pi/3 respectively, and
%                    must therefore be called without phi
%
% Outputs:
%
%     Q1,Q2,...    - quadrupolar coupling tensors for the two,
%                    three, four, or six sites as 3x3 matrices
%                    in Hz
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=weblab2nqi.m>

function varargout=weblab2nqi(C_q,eta_q,I,alpha,theta,phi)

% Four- and six-output modes take no phi
if nargin<6, phi=[]; end

% Check consistency
grumble(C_q,eta_q,I,alpha,theta,phi,nargin,nargout);

% Translate conventions and call eeqq2nqi
switch nargout
    
    case 2
        
        % Two output arguments
        varargout{1}=eeqq2nqi(C_q,eta_q,I,[-phi/2 theta alpha]);
        varargout{2}=eeqq2nqi(C_q,eta_q,I,[+phi/2 theta alpha]);
        
    case 3
        
        % Three output arguments
        varargout{1}=eeqq2nqi(C_q,eta_q,I,[-phi theta alpha]);
        varargout{2}=eeqq2nqi(C_q,eta_q,I,[ 0   theta alpha]);
        varargout{3}=eeqq2nqi(C_q,eta_q,I,[+phi theta alpha]);
        
    case 4
        
        % Four output arguments, fixed azimuth grid
        varargout{1}=eeqq2nqi(C_q,eta_q,I,[0*pi/2 theta alpha]);
        varargout{2}=eeqq2nqi(C_q,eta_q,I,[1*pi/2 theta alpha]);
        varargout{3}=eeqq2nqi(C_q,eta_q,I,[2*pi/2 theta alpha]);
        varargout{4}=eeqq2nqi(C_q,eta_q,I,[3*pi/2 theta alpha]);
        
    case 6
        
        % Six output arguments, fixed azimuth grid
        varargout{1}=eeqq2nqi(C_q,eta_q,I,[0*pi/3 theta alpha]);
        varargout{2}=eeqq2nqi(C_q,eta_q,I,[1*pi/3 theta alpha]);
        varargout{3}=eeqq2nqi(C_q,eta_q,I,[2*pi/3 theta alpha]);
        varargout{4}=eeqq2nqi(C_q,eta_q,I,[3*pi/3 theta alpha]);
        varargout{5}=eeqq2nqi(C_q,eta_q,I,[4*pi/3 theta alpha]);
        varargout{6}=eeqq2nqi(C_q,eta_q,I,[5*pi/3 theta alpha]);
        
    otherwise
        
        % Complain and bomb out
        error('incorrect number of output arguments');

end

end

% Consistency enforcement
function grumble(C_q,eta_q,I,alpha,theta,phi,n_ins,n_outs)
if ((n_outs==2)||(n_outs==3))&&(n_ins~=6)
    error('phi is required in the two- and three-output modes.');
end
if ((n_outs==4)||(n_outs==6))&&(n_ins~=5)
    error('phi is not accepted in the four- and six-output modes.');
end
if (~isnumeric(C_q))||(~isnumeric(eta_q))||(~isnumeric(I))||...
   (~isnumeric(alpha))||(~isnumeric(theta))||(~isnumeric(phi))
    error('all inputs must be numeric.');
end
if (~isreal(C_q))||(~isreal(eta_q))||(~isreal(I))||...
   (~isreal(alpha))||(~isreal(theta))||(~isreal(phi))
    error('all inputs must be real.');
end
if (~isscalar(C_q))||(~isscalar(eta_q))||(~isscalar(I))||...
   (~isscalar(alpha))||(~isscalar(theta))||...
   ((n_ins==6)&&(~isscalar(phi)))
    error('all inputs must be scalar.');
end
if (numel(I)~=1)||(I<1)||(mod(2*I+1,1)~=0)
    error('I must be an integer or half-integer greater or equal to 1.');
end
end

% Тот, кому отказывают в равенстве,
% вынужден добиваться превосходства.
%
% Российская Газета

