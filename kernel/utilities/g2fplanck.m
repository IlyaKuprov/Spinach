% Returns gradient operators within the Fokker-Planck formalism
% used in the imaging module of Spinach. Syntax:
%
%                G=g2fplanck(spin_system,parameters)
%
% Parameters:
%
%   parameters.dims    - a vector with one, two or three 
%                        elements giving the dimensions
%                        of the box, metres
%
%   parameters.npts    - a vector with one, two or three
%                        elements giving number of points
%                        in each dimension of the box
%
% Outputs:
%
%   G - a cell array with the three gradient operators 
%       ordered as {Gx,Gy,Gz}, normalised to 1 T/m, empty
%       matrices for non-exitent dimensions
%
% Note: gradients are assumed to be linear and centered on the
%       middle of the sample.
%
% Note: the direct product order is Z(x)Y(x)X(x)Spin, this cor-
%       responds to a column-wise vectorization of a 3D array
%       with dimensions ordered as [X Y Z].
%
% Note: polyadic objects are returned, use inflate() to get the
%       corresponding sparse matrix.
%
% a.j.allami@soton.ac.uk
% i.kuprov@soton.ac.uk
% m.g.concilio@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=g2fplanck.m>

function G=g2fplanck(spin_system,parameters)

% Check consistency
grumble(spin_system,parameters);

% Call Spinach to build magnet Zeeman Hamiltonian (per Tesla)
H=hamiltonian(assume(spin_system,'labframe','zeeman'))/spin_system.inter.magnet;

% Empty arrays as default
Gx=[]; Gy=[]; Gz=[]; %#ok<NASGU>

% Build gradient operators
switch numel(parameters.npts)
    
    case 1
        
        % Generate normalized gradient operators in 1D
        Gx=linspace(-0.5,0.5,parameters.npts(1));
        Gx=spdiags(Gx',0,parameters.npts(1),parameters.npts(1));
        Gx=parameters.dims(1)*polyadic({{Gx,H}});
        
    case 2
        
        % Generate normalized gradient operators in 2D
        Gx=parameters.dims(1)*linspace(-0.5,0.5,parameters.npts(1));
        Gy=parameters.dims(2)*linspace(-0.5,0.5,parameters.npts(2));
        Gx=spdiags(Gx',0,parameters.npts(1),parameters.npts(1));
        Gy=spdiags(Gy',0,parameters.npts(2),parameters.npts(2));
        Gx=polyadic({{opium(parameters.npts(2),1),Gx,H}}); 
        Gy=polyadic({{Gy,opium(parameters.npts(1),1),H}});
        
        % Rotate gradients if necessary
        if isfield(parameters,'grad_angles')
            R=[ cos(parameters.grad_angles) sin(parameters.grad_angles);
               -sin(parameters.grad_angles) cos(parameters.grad_angles)];
            Gx_new=Gx*R(1,1)+Gy*R(1,2);
            Gy_new=Gx*R(2,1)+Gy*R(2,2);
            Gx=Gx_new; Gy=Gy_new;
        end

    case 3
        
        % Generate normalized gradient operators
        Gx=parameters.dims(1)*linspace(-0.5,0.5,parameters.npts(1));
        Gy=parameters.dims(2)*linspace(-0.5,0.5,parameters.npts(2));
        Gz=parameters.dims(3)*linspace(-0.5,0.5,parameters.npts(3));
        Gx=spdiags(Gx',0,parameters.npts(1),parameters.npts(1));
        Gy=spdiags(Gy',0,parameters.npts(2),parameters.npts(2));
        Gz=spdiags(Gz',0,parameters.npts(3),parameters.npts(3));
        Gx=polyadic({{opium(parameters.npts(3),1),opium(parameters.npts(2),1),Gx,H}});
        Gy=polyadic({{opium(parameters.npts(3),1),Gy,opium(parameters.npts(1),1),H}});
        Gz=polyadic({{Gz,opium(parameters.npts(2),1),opium(parameters.npts(1),1),H}});
        
        % Rotate gradients if necessary
        if isfield(parameters,'grad_angles')
            R=euler2dcm(parameters.grad_angles);
            Gx_new=Gx*R(1,1)+Gy*R(1,2)+Gz*R(1,3);
            Gy_new=Gx*R(2,1)+Gy*R(2,2)+Gz*R(2,3);
            Gz_new=Gx*R(3,1)+Gy*R(3,2)+Gz*R(3,3);
            Gx=Gx_new; Gy=Gy_new; Gz=Gz_new;
        end
        
    otherwise
        
        % Complain and bomb out
        error('incorrect number of spatial dimensions.');
        
end

% Pack the operators
G={Gx,Gy,Gz};

end

% Consistency enforcement
function grumble(spin_system,parameters) %#ok<INUSD>
if spin_system.inter.magnet==0
    error('the primary magnet field must be non-zero.');
end
end

% There is nothing of any importance in life - except how well you do
% your work. Nothing. Only that. Whatever else you are, will come from
% that. It's the only measure of human value. All the codes of ethics
% they'll try to ram down your throat are just so much paper money put
% out by swindlers to fleece people of their virtues. The code of com-
% petence is the only system of morality that's on a gold standard. 
% When you grow up, you'll know what I mean.
%
% Ayn Rand, "Atlas Shrugged"

