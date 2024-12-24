% A basic hydrodynamics infrastructure provider, returns first
% derivative operators with respect to the three sample coordi-
% nates. Syntax:
%
%              [Fx,Fy,Fz]=hydrodynamics(parameters)
%
% Parameters:
%
%   parameters.dims    - dimensions of the sample (meters),
%                        one, two, or three-element row
%                        vector
%
%   parameters.npts    - number of points in each dimension
%                        of the sample, one, two, or three-
%                        element row vector
%
%   parameters.deriv   - {'fourier'} requests Fourier diffe-
%                        rentiation matrices; {'period',n}
%                        requests n-point central finite-
%                        difference matrices with periodic
%                        boundary conditions
%
% Outputs:
%
%   Fx, Fy, Fz         - derivative matrices, SI units
% 
% Note: the direct product order is Z(x)Y(x)X(x)Spin, this cor-
%       responds to a column-wise vectorization of a 3D array
%       with dimensions ordered as [X Y Z].
%
% Note: polyadic objects are returned, use inflate() to get the
%       corresponding sparse matrix.
%
% ilya.kuprov@weizmann.ac.uk
% a.j.allami@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=hydrodynamics.m>

function [Fx,Fy,Fz]=hydrodynamics(spin_system,parameters)

% Check consistency
grumble(parameters);

% Build derivative operators
switch parameters.deriv{1}
    
    case 'period'
        
        % Finite-difference derivatives
        if isscalar(parameters.npts)
            Dx=fdmat(parameters.npts(1),parameters.deriv{2},1)/(parameters.dims(1)/parameters.npts(1));
        end
        if numel(parameters.npts)==2
            Dx=fdmat(parameters.npts(1),parameters.deriv{2},1)/(parameters.dims(1)/parameters.npts(1));
            Dy=fdmat(parameters.npts(2),parameters.deriv{2},1)/(parameters.dims(2)/parameters.npts(2));
        end
        if numel(parameters.npts)==3
            Dx=fdmat(parameters.npts(1),parameters.deriv{2},1)/(parameters.dims(1)/parameters.npts(1));
            Dy=fdmat(parameters.npts(2),parameters.deriv{2},1)/(parameters.dims(2)/parameters.npts(2));
            Dz=fdmat(parameters.npts(3),parameters.deriv{2},1)/(parameters.dims(3)/parameters.npts(3));
        end
        
   case 'fourier'
        
        % Fourier derivatives
        if isscalar(parameters.npts)
            [~,Dx]=fourdif(parameters.npts(1),1); Dx=(2*pi/parameters.dims(1))*Dx;
        end
        if numel(parameters.npts)==2
            [~,Dx]=fourdif(parameters.npts(1),1); Dx=(2*pi/parameters.dims(1))*Dx;
            [~,Dy]=fourdif(parameters.npts(2),1); Dy=(2*pi/parameters.dims(2))*Dy;
        end
        if numel(parameters.npts)==3
            [~,Dx]=fourdif(parameters.npts(1),1); Dx=(2*pi/parameters.dims(1))*Dx;
            [~,Dy]=fourdif(parameters.npts(2),1); Dy=(2*pi/parameters.dims(2))*Dy;
            [~,Dz]=fourdif(parameters.npts(3),1); Dz=(2*pi/parameters.dims(3))*Dz;
        end
        
    otherwise
        
        % Complain and bomb out
        error('unrecognized derivative operator type.');
        
end

% Kron up derivative operators
Fx=[]; Fy=[]; Fz=[];
if isscalar(parameters.npts)
    Fx=-1i*polyadic({{Dx}}); 
end
if numel(parameters.npts)==2
    Fx=-1i*polyadic({{opium(parameters.npts(2),1),Dx}});
    Fy=-1i*polyadic({{Dy,opium(parameters.npts(1),1)}});
end
if numel(parameters.npts)==3
    Fx=-1i*polyadic({{opium(parameters.npts(3),1),opium(parameters.npts(2),1),Dx}});
    Fy=-1i*polyadic({{opium(parameters.npts(3),1),Dy,opium(parameters.npts(1),1)}});
    Fz=-1i*polyadic({{Dz,opium(parameters.npts(2),1),opium(parameters.npts(1),1)}});
end

% Inflate derivative operators
if ~ismember('polyadic',spin_system.sys.enable)
   Fx=inflate(Fx); Fy=inflate(Fy); Fz=inflate(Fz);
end

end

% Consistency enforcement
function grumble(parameters)
if ~isfield(parameters,'dims')
    error('sample dimensions must be specified in parameters.dims variable.');
end
if (~isnumeric(parameters.dims))||(~isreal(parameters.dims))||...
   (any(parameters.dims<=0))||(numel(parameters.dims)<1)||...
   (numel(parameters.dims)>3)||(size(parameters.dims,1)~=1)
    error('parameters.dims must be a 1, 2, or 3-element row vector of positive real numbers.');
end
if ~isfield(parameters,'npts')
    error('grid sizes must be specified in parameters.npts variable.');
end
if (~isnumeric(parameters.npts))||(~isreal(parameters.npts))||...
   (numel(parameters.npts)<1)||(numel(parameters.npts)>3)||...
   (size(parameters.npts,1)~=1)||any(mod(parameters.npts,1))
    error('parameters.dims must be a 1, 2, or 3-element row vector of positive real numbers.');
end
if numel(parameters.npts)~=numel(parameters.dims)
    error('parameters.dims and parameters.npts must have the same number of elements.');
end
if ~isfield(parameters,'deriv')
    error('differentiation method must be specified in parameters.deriv variable.');
end
if (~iscell(parameters.deriv))||(numel(parameters.deriv)<1)||(numel(parameters.deriv)>2)
    error('parameters.deriv must be a cell array with one or two elements');
end
if (~ischar(parameters.deriv{1}))||(~ismember(parameters.deriv{1},{'period','fourier'}))
    error('the first element of parameters.deriv must be ''period'' or ''fourier''.');
end
if strcmp(parameters.deriv{1},'period')
    if (~isnumeric(parameters.deriv{2}))||(~isreal(parameters.deriv{2}))||...
       (numel(parameters.deriv{2})~=1)||mod(parameters.deriv{2},1)
        error('stencil size in the second element of parameters.deriv must be a positive integer.');
    end
    if parameters.deriv{2}>7
        error('differentiation stencil size greater than 7 is not a good idea - use a bigger grid.');
    end
end
if strcmp(parameters.deriv{1},'fourier')
    if numel(parameters.deriv)>1
        error('stencil size parameter does not apply to Fourier differentiation.');
    end
end
if any(parameters.npts<10)
    error('a spatial grid with fewer than 10 points is not a good idea - use a bigger grid.');
end
end   

% The weather vane was nailed down, and the wind was 
% blowing forlornly in the direction indicated.
%
% Viktor Shenderovich

