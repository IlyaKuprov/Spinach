% Converts the Weblab one-cone model parameters (see weblab_cone.png)
% into NQI tensors used by Spinach. Syntax:
%
%         [Q1,Q2,...]=weblab2nqi(C_q,eta_q,I,alpha,theta,phi)
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
%                    (see weblab_cone.png), in radians
%
% Outputs:
%
%     Q1,Q2,...    - quadrupolar coupling tensors for the two,
%                    three, four, or six sites as 3x3 matrices
%                    in Hz
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=weblab2nqi.m>

function varargout=weblab2nqi(C_q,eta_q,I,alpha,theta,phi)

% Check consistency
grumble(C_q,eta_q,I,alpha,theta,phi);

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
        
        % Four output arguments, fixed phi values
        if ~isempty(phi)
            disp('WARNING: phi was overwritten to [0 1 2 3]*pi/4');
        end
        varargout{1}=eeqq2nqi(C_q,eta_q,I,[0*pi/2 theta alpha]);
        varargout{2}=eeqq2nqi(C_q,eta_q,I,[1*pi/2 theta alpha]);
        varargout{3}=eeqq2nqi(C_q,eta_q,I,[2*pi/2 theta alpha]);
        varargout{4}=eeqq2nqi(C_q,eta_q,I,[3*pi/2 theta alpha]);
        
    case 6
        
        % Six output arguments, fixed phi values
        if ~isempty(phi)
            disp('WARNING: phi was overwritten to [0 1 2 3 4 5]*pi/6');
        end
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
function grumble(C_q,eta_q,I,alpha,theta,phi)
if (~isnumeric(C_q))||(~isnumeric(eta_q))||(~isnumeric(I))||...
   (~isnumeric(alpha))||(~isnumeric(theta))||(~isnumeric(phi))
    error('all inputs must be numeric.');
end
if (~isreal(C_q))||(~isreal(eta_q))||(~isreal(I))||...
   (~isreal(alpha))||(~isreal(theta))||(~isreal(phi))
    error('all inputs must be real.');
end
if (~isscalar(C_q))||(~isscalar(eta_q))||(~isscalar(I))||...
   (~isscalar(alpha))||(~isscalar(theta))||(~isscalar(phi))
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

