% Translates a stationary 3D velocity field and a diffusion tensor
% field into a Fokker-Planck evolution generator. Syntax:
%
%                F=v2fplanck(spin_system,parameters)
%
% The following parameters are necessary:
%
%   parameters.u       - X components of the velocity vectors
%                        for each voxel in the sample, m/s
%
%   parameters.v       - Y components of the velocity vectors
%                        for each voxel in the sample, m/s
%
%   parameters.w       - Z components of the velocity vectors
%                        for each voxel in the sample, m/s
%
%   parameters.diff    - diffusion coefficient or 3x3 tensor, m^2/s
%                        for situations when this parameter is the 
%                        same in every voxel
%
%   parameters.dxx     - Cartesian components of the diffusion
%   parameters.dxy       tensor for each voxel of the sample
%        ...
%   parameters.dzz
%
%   parameters.dims    - dimensions of the 3D box, meters
%
%   parameters.npts    - number of points in each dimension
%                        of the 3D box
%
%   parameters.deriv   - {'fourier'} uses Fourier diffe-
%                        rentiation matrices; {'period',n}
%                        requests n-point central finite-
%                        difference matrices with periodic
%                        boundary conditions
%
% Outputs:
%
%   F - spatial dynamics generator
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
%
% <https://spindynamics.org/wiki/index.php?title=v2fplanck.m>

function F=v2fplanck(spin_system,parameters)

% Check consistency
grumble(parameters);

% Get translation generators
[Fx,Fy,Fz]=hydrodynamics(spin_system,parameters);

% Count voxels
nvoxels=prod(parameters.npts);

% Default is no flow on X
if ~isfield(parameters,'u'), parameters.u=0; end

% Flow in X direction
if isscalar(parameters.u)
    
    % Uniform flow on X
    F=parameters.u*Fx;
    
else
    
    % Non-uniform flow on X
    F=spdiags(Fx*parameters.u(:),0,nvoxels,nvoxels)+...
      spdiags(parameters.u(:),0,nvoxels,nvoxels)*Fx;
    
end

% Uniform X diffusion
if isfield(parameters,'diff')
    
    % Diffusion type
    if isscalar(parameters.diff)
    
        % Add isotropic diffusion on X
        F=F-1i*parameters.diff*Fx*Fx;
        
    else
        
        % Add anisotropic diffusion on X
        F=F-1i*parameters.diff(1,1)*Fx*Fx;
        
    end
    
end

% Non-uniform X diffusion
if isfield(parameters,'dxx')
    
    % Add non-uniform diffusion on X
    F=F-1i*Fx*spdiags(parameters.dxx(:),0,nvoxels,nvoxels)*Fx;
    
end

% When sample has Y dimension
if numel(parameters.npts)>1
    
    % Default is no flow on Y
    if ~isfield(parameters,'v'), parameters.v=0; end

    % Flow in Y direction
    if isscalar(parameters.v)
        
        % Add uniform flow on Y
        F=F+parameters.v*Fy;
        
    else
        
        % Add non-uniform flow on Y
        F=F+spdiags(Fy*parameters.v(:),0,nvoxels,nvoxels)+...
            spdiags(parameters.v(:),0,nvoxels,nvoxels)*Fy;
        
    end
    
    % Uniform diffusion in XY
    if isfield(parameters,'diff')
        
        % Diffusion type
        if isscalar(parameters.diff)
    
            % Add isotropic diffusion on Y
            F=F-1i*parameters.diff*Fy*Fy;
            
        else
            
            % Update anisotropic diffusion to XY
            F=F-1i*parameters.diff(1,2)*Fx*Fy;
            F=F-1i*parameters.diff(2,1)*Fy*Fx;
            F=F-1i*parameters.diff(2,2)*Fy*Fy;
            
        end
        
    end
    
    % Non-uniform diffusion in XY 
    if all(isfield(parameters,{'dxx','dxy','dyx','dyy'}))
    
        % Update non-uniform diffusion to XY
        F=F-1i*(Fx*spdiags(parameters.dxy(:),0,nvoxels,nvoxels)*Fy+...
                Fy*spdiags(parameters.dyx(:),0,nvoxels,nvoxels)*Fx+...
                Fy*spdiags(parameters.dyy(:),0,nvoxels,nvoxels)*Fy);
    
    end
    
end

% When sample has Z dimension
if numel(parameters.npts)>2
    
    % Default is no flow on Z
    if ~isfield(parameters,'w'), parameters.w=0; end

    % Flow in Z direction
    if isscalar(parameters.w)
        
        % Add uniform flow on Z
        F=F+parameters.w*Fz;
        
    else
        
        % Add non-uniform flow on Z
        F=F+spdiags(Fz*parameters.w(:),0,nvoxels,nvoxels)+...
            spdiags(parameters.w(:),0,nvoxels,nvoxels)*Fz;
        
    end
    
    % Uniform diffusion in XYZ
    if isfield(parameters,'diff')
        
        % Diffusion type
        if isscalar(parameters.diff)
    
            % Add isotropic diffusion on Z
            F=F-1i*parameters.diff*Fz*Fz;
            
        else
            
            % Update anisotropic diffusion to XYZ
            F=F-1i*parameters.diff(1,3)*Fx*Fz;
            F=F-1i*parameters.diff(2,3)*Fy*Fz;
            F=F-1i*parameters.diff(3,3)*Fz*Fz;
            F=F-1i*parameters.diff(3,2)*Fz*Fy;
            F=F-1i*parameters.diff(3,1)*Fz*Fx;
            
        end
        
    end
            
    % Non-uniform diffusion in XYZ
    if all(isfield(parameters,{'dxx','dxy','dxz',...
                               'dyx','dyy','dyz',...
                               'dzx','dzy','dzz'}))
    
        % Update non-uniform diffusion to XYZ
        F=F-1i*(Fx*spdiags(parameters.dxz(:),0,nvoxels,nvoxels)*Fz+...
                Fy*spdiags(parameters.dyz(:),0,nvoxels,nvoxels)*Fz+...
                Fz*spdiags(parameters.dzx(:),0,nvoxels,nvoxels)*Fx+...
                Fz*spdiags(parameters.dzy(:),0,nvoxels,nvoxels)*Fy+...
                Fz*spdiags(parameters.dzz(:),0,nvoxels,nvoxels)*Fz);
    
    end
    
end

% Clean up the array
F=clean_up(spin_system,F,spin_system.tols.liouv_zero);

% Kron up with the spin
spn_dim=size(spin_system.bas.basis,1);
F=kron(F,opium(spn_dim,1));

end

% Consistency enforcement
function grumble(parameters)
if ~isfield(parameters,'npts')
    error('the number of points in the spatial grid must be specified in parameters.npts field.');
end
if (~isnumeric(parameters.npts))||(~isreal(parameters.npts))||...
   (~isrow(parameters.npts))||any(parameters.npts<1)||...
   (numel(parameters.npts)>3)||any(mod(parameters.npts,1)~=0)
    error('parameters.npts must be a 1, 2 or 3-element row vector of positive integers.');
end
if isfield(parameters,'diff')&&any(isfield(parameters,...
   {'dxx','dxy','dxz','dyx','dyy','dyz','dzx','dzy','dzz'}))
    error('specify either parameters.diff, or parameters.d**, but not both.');
end
if isscalar(parameters.npts)
    if isfield(parameters,'v')||isfield(parameters,'w')
        error('parameters.v,w do not apply to one-dimensional samples.');
    end
    if any(isfield(parameters,{'dxy','dxz','dyx','dyy','dyz','dzx','dzy','dzz'}))
        error('parameters.dxy,dxz,dyx,dyy,dyz,dzx,dzy,dzz do not apply to one-dimensional samples.');
    end
    if isfield(parameters,'u')
        if (~isnumeric(parameters.u))||(~isreal(parameters.u))||...
           (~iscolumn(parameters.u))||any(~isfinite(parameters.u))
            error('parameters.u must be a real column vector.');
        end
        if numel(parameters.u)~=parameters.npts
            error('the number of elements in parameters.u must be equal to parameters.npts');
        end
        for n={'dxx'}
            if isfield(parameters,n{1})
                d=getfield(parameters,n{1}); %#ok<GFLD>
                if (~isnumeric(d))||(~isreal(d))||...
                   (~iscolumn(d))||any(~isfinite(d))
                    error(['parameters.' n{1} ' must be a real column vector.']);
                end
            end
        end
    end
end
if numel(parameters.npts)==2
    if isfield(parameters,'w')
        error('parameters.w does not apply to two-dimensional samples.');
    end
    if any(isfield(parameters,{'dxz','dyz','dzx','dzy','dzz'}))
        error('parameters.dxz,dyz,dzx,dzy,dzz do not apply to two-dimensional samples.');
    end
    if isfield(parameters,'u')
        if (~isnumeric(parameters.u))||(~isreal(parameters.u))||...
           any(size(parameters.u)~=parameters.npts)||any(~isfinite(parameters.u(:)))
            error(['parameters.u must be a real array of dimension ' num2str(parameters.npts)]);
        end
    end
    if isfield(parameters,'v')
        if (~isnumeric(parameters.v))||(~isreal(parameters.v))||...
           any(size(parameters.v)~=parameters.npts)||any(~isfinite(parameters.v(:)))
            error(['parameters.v must be a real array of dimension ' num2str(parameters.npts)]);
        end
    end
    for n={'dxx','dxy','dyx','dyy'}
        if isfield(parameters,n{1})
            d=getfield(parameters,n{1}); %#ok<GFLD>
            if (~isnumeric(d))||(~isreal(d))||any(size(d)~=parameters.npts)||...
                any(~isfinite(d(:)))||any(d(:)<0)
                error(['parameters.' n{1} ' must be a real array of dimension ' num2str(parameters.npts)]);
            end
        end
    end
    if any(isfield(parameters,{'dxx','dxy','dyx','dyy'}))&&...
       (~all(isfield(parameters,{'dxx','dxy','dyx','dyy'})))
        error('parameters.dxx,dxy,dyx,dyy must be specified simultaneously.');
    end
end
if (numel(parameters.npts)==3)
    if isfield(parameters,'u')
        if (~isnumeric(parameters.u))||(~isreal(parameters.u))||...
           any(size(parameters.u)~=parameters.npts)||any(~isfinite(parameters.u(:)))
            error(['parameters.u must be a real array of dimension ' num2str(parameters.npts)]);
        end
    end
    if isfield(parameters,'v')
        if (~isnumeric(parameters.v))||(~isreal(parameters.v))||...
           any(size(parameters.v)~=parameters.npts)||any(~isfinite(parameters.v(:)))
            error(['parameters.v must be a real array of dimension ' num2str(parameters.npts)]);
        end
    end
    if isfield(parameters,'w')
        if (~isnumeric(parameters.w))||(~isreal(parameters.w))||...
           any(size(parameters.w)~=parameters.npts)||any(~isfinite(parameters.w(:)))
            error(['parameters.w must be a real array of dimension ' num2str(parameters.npts)]);
        end
    end
    for n={'dxx','dxy','dxz','dyx','dyy','dyz','dzx','dzy','dzz'}
        if isfield(parameters,n{1})
            d=getfield(parameters,n{1}); %#ok<GFLD>
            if (~isnumeric(d))||(~isreal(d))||any(size(d)~=parameters.npts)||...
                any(~isfinite(d(:)))||any(d(:)<0)
                error(['parameters.' n{1} ' must be a real array of dimension ' num2str(parameters.npts)]);
            end
        end
    end
    if any(isfield(parameters,{'dxx','dxy','dxz','dyx','dyy','dyz','dzx','dzy','dzz'}))&&...
       (~all(isfield(parameters,{'dxx','dxy','dxz','dyx','dyy','dyz','dzx','dzy','dzz'})))
        error('parameters.dxx,dxy,dxz,dyx,dyy,dyz,dzx,dzy,dzz must be specified simultaneously.');
    end
end
if isfield(parameters,'diff')
    if (~isnumeric(parameters.diff))||(~isreal(parameters.diff))
        error('parameters.diff must a real scalar or matrix.');
    end
end
if ~isfield(parameters,'dims')
    error('sample dimensions must be specified in parameters.dims field.');
end
if (~isnumeric(parameters.dims))||(~isreal(parameters.dims))||...
   (any(~isfinite(parameters.dims)))||any(parameters.dims<=0)
    error('parameters.dims must be a row vector of positive real numbers.');
end
if numel(parameters.dims)~=numel(parameters.npts)
    error('the number of elements in parameters.dims and parameters.npts must be the same.');
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

% "He's spending a year dead for tax reasons."
%
% Douglas Adams, The Hitchhiker's Guide to the Galaxy

