% Draws shielding tensors and their eigensystems. Two styles are
% implemented:
%
%   A. Ellipsoids (symmetric tensors only):
%
%      1. A unit sphere in a Cartesian space is scaled by 
%         abs(Axx) in the x direction, abs(Ayy) in the y 
%         direction and abs(Azz) in the z direction, where
%         Axx, Ayy, Azz are the eigenvalues of the CST ten-
%         sor in units of ppm.
%
%      2. A set of axes is drawn inside the sphere with a 
%         red axis for a positive eigenvalue, and a blue
%         axis for a negative one.
%
%      3. The sphere is translated to the point of corres-
%         ponding nucleus and rotated into the molecular 
%         frame of reference.
%
%   B. Spherical harmonics (default):
%
%      1. The matrix is converted into irreducible spheri-
%         cal tensor operator coefficients.
%
%      2. The coefficients are placed in front of the cor-
%         responding spherical harmonics, which are plot-
%         ted in three dimensions and translated to the po-
%         int of the corresponding nucleus.
%
% Syntax:
%
%       cst_display(props,atoms,scaling,conmatrix,options)
%
% Arguments:
%
%               props - output of gparse function
%
%               atoms - a cell array of element symbols
%                       or a vector of integers, indica-
%                       ting the atoms for which shiel-
%                       ding tensors should be visuali-
%                       sed, e.g. {'C','H'} or [1 2 5]
%
%             scaling - a factor to scale the tensors
%                       by for visualisation
%
%           conmatrix - binary connectivity matrix, 1
%                       if a pair of atoms should be
%                       connected by a bond. If an em-
%                       pty vector is supplied, 1.6 
%                       Angstrom cutoff distance is used
%
%       options.style - 'ellipsoids' or 'harmonics'
%
%    options.kill_iso - set to true() to eliminate the
%                       isotropic parts of tensors be-
%                       fore plotting
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=cst_display.m>

function cst_display(props,atoms,scaling,conmatrix,options)

% Check consistency
grumble(props,atoms,scaling,conmatrix);

% Set defaults
if (~exist('options','var'))||...
   (~isfield(options,'style'))
   options.style='harmonics';
end
if (~exist('options','var'))||...
   (~isfield(options,'kill_iso'))
   options.kill_iso=false();
end
    
% Set up graphics
light('Position',[-2, 2,20]); 
light('Position',[10,10,10]); hold('on');

% Draw the molecule
molplot(props.std_geom,conmatrix);

% Sample a unit sphere
k=5; npts=2^k-1; 
theta=linspace(0,pi,npts);
phi=linspace(0,2*pi,npts);
[T,P]=ndgrid(theta,phi);
X=sin(T).*cos(P); X=X(:)';
Y=sin(T).*sin(P); Y=Y(:)';
Z=cos(T);         Z=Z(:)';

% Atoms by number or name
if isnumeric(atoms)
    idx=atoms(:)';
elseif iscell(atoms)
    idx=find(ismember(props.symbols,atoms(:)));
else
    error('incorrect atoms specification.');
end

% Loop over atoms
for n=idx

    % Get the current tensor
    cst=props.cst{n};
        
    % Eliminate the isotroic part
    if options.kill_iso
        cst=cst-eye(3)*trace(cst)/3;
    end
        
    % Pick the style
    switch options.style
          
        % Ellipsoid plots
        case 'ellipsoids'
                
            % Diagonalize the tensor
            [eigvecs,eigvals]=eig(cst,'vector');
                
            % Make sure the tensor is symmetric
            if norm(eigvecs'*eigvecs-eye(3),2)>1e-3
                error('significant antisymmetry found, use options.style=''harmonics''');
            end
            
            % Scale vertex coordinates
            coords=[X*eigvals(1)*scaling;
                    Y*eigvals(2)*scaling;
                    Z*eigvals(3)*scaling];
            
            % Rotate vertex coordinates
            coords=eigvecs*coords;

            % Translate vertex coordinates
            coords=coords+props.std_geom(n,:)';

            % Reshape for plotting
            Xp=reshape(coords(1,:),[npts npts]);
            Yp=reshape(coords(2,:),[npts npts]);
            Zp=reshape(coords(3,:),[npts npts]);

            % Simple grey colour
            Cp=ones([npts npts 3])/2;

            % Draw the surface
            surf(Xp,Yp,Zp,Cp,'FaceAlpha',0.5,'LineStyle','none',...
                             'FaceColor','flat','FaceLighting','flat');

            % Scale eigenvectors
            vector_a=eigvecs(:,1)*eigvals(1)*scaling;
            vector_b=eigvecs(:,2)*eigvals(2)*scaling;
            vector_c=eigvecs(:,3)*eigvals(3)*scaling;

            % Color eigenvectors
            col_a='r-'; col_b='r-'; col_c='r-';
            if eigvals(1)<0, col_a='b-'; end
            if eigvals(2)<0, col_b='b-'; end
            if eigvals(3)<0, col_c='b-'; end

            % Draw eigenvectors
            plot3([props.std_geom(n,1)-vector_a(1) props.std_geom(n,1)+vector_a(1)],...
                  [props.std_geom(n,2)-vector_a(2) props.std_geom(n,2)+vector_a(2)],...
                  [props.std_geom(n,3)-vector_a(3) props.std_geom(n,3)+vector_a(3)],col_a,'LineWidth',1);
            plot3([props.std_geom(n,1)-vector_b(1) props.std_geom(n,1)+vector_b(1)],...
                  [props.std_geom(n,2)-vector_b(2) props.std_geom(n,2)+vector_b(2)],...
                  [props.std_geom(n,3)-vector_b(3) props.std_geom(n,3)+vector_b(3)],col_b,'LineWidth',1);
            plot3([props.std_geom(n,1)-vector_c(1) props.std_geom(n,1)+vector_c(1)],...
                  [props.std_geom(n,2)-vector_c(2) props.std_geom(n,2)+vector_c(2)],...
                  [props.std_geom(n,3)-vector_c(3) props.std_geom(n,3)+vector_c(3)],col_c,'LineWidth',1);

            % Spherical harmonic plots
        case 'harmonics'

            % Get IST coefficients
            [rank0,rank1,rank2]=mat2sphten(cst);

            % Compute radii
            R=real(rank0   *spher_harmon(0,0,T,P) +...
                   rank1(1)*spher_harmon(1,+1,T,P)+...
                   rank1(2)*spher_harmon(1, 0,T,P)+...
                   rank1(3)*spher_harmon(1,-1,T,P)+...
                   rank2(1)*spher_harmon(2,+2,T,P)+...
                   rank2(2)*spher_harmon(2,+1,T,P)+...
                   rank2(3)*spher_harmon(2, 0,T,P)+...
                   rank2(4)*spher_harmon(2,-1,T,P)+...
                   rank2(5)*spher_harmon(2,-2,T,P)); R=R(:)';

            % Scale vertex coordinates
            coords=scaling*[R.*X; R.*Y; R.*Z];

            % Translate vertex coordinates
            coords=coords+props.std_geom(n,:)';

            % Reshape for plotting
            Xp=reshape(coords(1,:),[npts npts]);
            Yp=reshape(coords(2,:),[npts npts]);
            Zp=reshape(coords(3,:),[npts npts]);

            % Red for pos, blue for neg
            Cp=zeros([npts npts 3]);
            Cp(:,:,1)=(reshape(R,[npts npts])>0)/2;
            Cp(:,:,3)=(reshape(R,[npts npts])<0)/2;

            % Draw the surface
            surf(Xp,Yp,Zp,Cp,'FaceAlpha',0.25,'LineStyle','none',...
                             'FaceColor','flat','FaceLighting','flat');

        otherwise

            % Complain and bomb out
            error('unknown plotting style');

    end

end

% Tidy up the picture
axis square; axis tight; axis equal;
set(gca,'Projection','perspective',...
        'XTickLabel',[],'YTickLabel',[],...
        'ZTickLabel',[],'TickLength',[0 0]); 
box on; kgrid;

end

% Consistency enforcement
function grumble(props,atoms,scaling,conmatrix)
if ~isfield(props,'std_geom')
    error('Gaussian data structure does not contain std_geom field.');
end
if ~isfield(props,'symbols')
    error('Gaussian data structure does not contain symbols field.');
end
if ~isfield(props,'cst')
    error('props structure does not contain cst field with chemical shielding tensors.');
end
if iscell(atoms)
    if ~all(cellfun(@ischar,atoms),'all')
        error('when a cell array, atoms must contain character strings.');
    end
elseif isnumeric(atoms)
    if (~isreal(atoms))||any(atoms<1,'all')||any(mod(atoms,1),'all')
        error('when a vector, atoms must contain positive integers.');
    end
    if any(atoms>numel(props.symbols),'all')
        error('specified index exceeds the number of atoms.');
    end
else
    error('atoms must be a vector of integers or a cell array of strings.');
end
if (~isnumeric(scaling))||(~isreal(scaling))||...
   (~isscalar(scaling))||(scaling<=0)
    error('scaling must be a positive real number.');
end
if ~isempty(conmatrix)
    if (size(conmatrix,1)~=size(conmatrix,2))||...
       (size(conmatrix,1)~=size(props.std_geom,1))
        error(['conmatrix must be a logical square matrix '...
               'of the same dimension as the number of atoms.']);
    end
end
end

% I'll be more enthusiastic about encouraging 
% thinking outside the box when there's eviden-
% ce of any thinking going on inside it.
%
% Terry Pratchett

