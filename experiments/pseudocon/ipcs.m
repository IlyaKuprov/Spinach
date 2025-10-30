% Solves the inverse problem for pseudocontact shift by recovering the
% source term in the Kuprov equation using Tikhonov regularization 
% procedure. Syntax: 
%
%     [source_cube,ranges,pred_pcs,err_ls,reg_a,reg_b]=...
%                   ipcs(parameters,nxyz,expt_pcs,chi,npoints,...
%                        lambda,margins,box_centre,box_size)
%
% Parameters:
% 
%     nxyz     - nuclear coordinates as [x y z] with multiple rows 
%                at which PCS has been measured, in Angstroms.
%
%     expt_pcs - pseudocontact shift in ppm at each nucleus.
%
%     chi      - electron magnetic susceptibility tensor, in units
%                of Angstrom^3.
%
%     npoints  - number of points in each dimension of the source
%                cube, a positive integer greater than 10.
%
%     lambda   - regularization parameters, the first element is 
%                the coefficient in front of the contrast term
%                and the second element is the coefficient in
%                front of the Tikhonov term.
%
%     margins  - a six-element vector specifying margins to take
%                around the bounding box of the nuclear coordina-
%                tes supplied, to account for the possibility that
%                the electron may be located on the periphery.
%
%     box_centre - a three-element vector in Angstrom specifying 
%                  the centre of the solution box
%
%     box_size   - a three-element vector in Angstrom specifying
%                  the size of the solution box
%
%     equation   - 'poisson' to recover the right hand side of the
%                   Poisson's equation, 'kuprov' to recover the
%                   probability density.
%
%     gpu        - set to 1 to enable GPU processing (much faster)
%
% Outputs:
% 
%     source_cube - source term cube with dimensions ordered as
%                   [X Y Z]
%
%     ranges      - Cartesian axis extents for the source cube as
%                   [xmin xmax ymin ymax zmin zmax] in Angstroms
%
%     pred_pcs    - pseudocontact shifts produced by the source
%                   cube returned in the first parameter
%
%     ls_err      - least squares error in ppm^2
%
%     reg_a       - contrast penalty term
%
%     reg_b       - Tikhonov penalty term the second element is the entropy
%                   penalty in the error functional, the third
%                   element is the tikhonov penalty in the error
%                   functional
%
% Note: for further information on the equations and algorithms used
%       in this function see http://dx.doi.org/10.1039/C4CP03106G
%
% gareth.charnock@oerc.ox.ac.uk
% ilya.kuprov@weizmann.ac.il
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=Ipcs.m>

function [source_cube,ranges,pred_pcs,err_ls,reg_a,reg_b]=ipcs(parameters,npoints,lambda)

% Check consistency
grumble(parameters,npoints,lambda);

% Allow GPUs Matlab has not yet seen
parallel.gpu.enableCUDAForwardCompatibility(true);

% Validate GPU option
if parameters.gpu
    parameters.gpu=logical(gpuDeviceCount); 
end

% Conserve GPU memory
if parameters.gpu
    G=gpuDevice(); reset(G);
    G.CachePolicy='minimum';
end

% Set internal scaling parameters
switch parameters.equation
    
    case 'kuprov'
        
        % Pin the density roughly around [0,1]
        probdens_scaling=1.0e3;
        tikhonov_scaling=0.640/npoints;
        contrast_scaling=1.0e3/npoints^3;
        
    case 'poisson'
        
        % Pin the density roughly around [0,1]
        probdens_scaling=1.0e4;
        tikhonov_scaling=0.640/npoints;
        contrast_scaling=1.0e3/npoints^3;
        
end

% Set the GPU data types
if parameters.gpu
    data_type='gpuArray';
else
    data_type='double';
end

% Choose appropriate axis ranges
ranges=[min(parameters.xyz(:,1))-parameters.margins(1) ...
        max(parameters.xyz(:,1))+parameters.margins(2) ...
        min(parameters.xyz(:,2))-parameters.margins(3) ...
        max(parameters.xyz(:,2))+parameters.margins(4) ...
        min(parameters.xyz(:,3))-parameters.margins(5) ...
        max(parameters.xyz(:,3))+parameters.margins(6)];
   
% Report to the user
disp(['X axis extents: ' num2str(ranges([1 2])) ' Angstrom.']);
disp(['Y axis extents: ' num2str(ranges([3 4])) ' Angstrom.']);
disp(['Z axis extents: ' num2str(ranges([5 6])) ' Angstrom.']);
disp(['Cube dimension: ' num2str([npoints npoints npoints]) ' points.']);

% Get the nuclear point sampling matrix
disp('Building nuclear sampling matrix...');
P=interpmat([npoints npoints npoints],ranges,parameters.xyz);

% Get the coordinate grid
[X,Y,Z]=ndgrid(linspace(ranges(1),ranges(2),npoints),...
               linspace(ranges(3),ranges(4),npoints),...
               linspace(ranges(5),ranges(6),npoints));

% Get the solution box and apply distance cutoffs
soln_box=(X>(parameters.box_cent(1)-parameters.box_size(1)/2))&...
         (X<(parameters.box_cent(1)+parameters.box_size(1)/2))&...
         (Y>(parameters.box_cent(2)-parameters.box_size(2)/2))&...
         (Y<(parameters.box_cent(2)+parameters.box_size(2)/2))&...
         (Z>(parameters.box_cent(3)-parameters.box_size(3)/2))&...
         (Z<(parameters.box_cent(3)+parameters.box_size(3)/2));
     
% Confine the solution to the relevant layer
if isfield(parameters,'confine')
    disp('Building the variational volume...');
    reject_region=false(size(soln_box));
    accept_region=false(size(soln_box));
    parfor m=1:size(parameters.xyz_all,1)
        distance=sqrt((X-parameters.xyz_all(m,1)).^2+...
                      (Y-parameters.xyz_all(m,2)).^2+...
                      (Z-parameters.xyz_all(m,3)).^2); %#ok<PFBNS>
        reject_region=reject_region|(distance<parameters.confine(1));
        accept_region=accept_region|(distance<parameters.confine(2));
    end
    soln_box=soln_box&(~reject_region)&accept_region;
end

% Pre-compute atom connectivity matrix
if isfield(parameters,'xyz_all')
    conmatrix=conmat(parameters.xyz_all,1.6);
end

% Index the active space
active_point_index=soln_box(:);

% Report to the user
disp([num2str(nnz(active_point_index)) ' points in the probability density']);

% Get the initial guess
if isfield(parameters,'guess')
    
    % If a guess is provided, resample onto the current grid
    disp('Interpolating the guess onto the current grid...');
    [X0,Y0,Z0]=ndgrid(linspace(ranges(1),ranges(2),size(parameters.guess,1)),...
                      linspace(ranges(3),ranges(4),size(parameters.guess,2)),...
                      linspace(ranges(5),ranges(6),size(parameters.guess,3)));
    guess=probdens_scaling*interpn(X0,Y0,Z0,parameters.guess,X,Y,Z,'cubic');
    
    % Eliminate resampling artefacts
    if strcmp(parameters.equation,'kuprov'), guess(guess<0)=0; end
    
    % Extract relevant points
    guess=guess(active_point_index);
    
else
    
    % If a guess is not provided, create one
    switch parameters.equation
        
        case 'kuprov'
            
            % Exact normalisation multiplier
            normint=(ranges(2)-ranges(1))*(ranges(4)-ranges(3))*(ranges(6)-ranges(5))*...
                     nnz(active_point_index)/probdens_scaling/npoints^3;
            
            % Uniform positive initial guess for Kuprov equation
            disp('Uniform initial guess...');
            guess=ones(nnz(active_point_index),1)/normint;
            
        case 'poisson'
            
            % Approximate normalisation multiplier
            normint=(ranges(2)-ranges(1))*(ranges(4)-ranges(3))*(ranges(6)-ranges(5))*...
                     nnz(active_point_index)/probdens_scaling/npoints^3;
            
            % Random initial guess for Poisson equation
            disp('Random initial guess...');
            guess=rand(nnz(active_point_index),1)/normint;
            
    end
       
end

% Generate Fourier derivative multipliers
[X,Y,Z]=ndgrid(fftdiff(1,npoints,1),fftdiff(1,npoints,1),fftdiff(1,npoints,1));

% Get Laplace operator in Fourier space
L=X.^2+Y.^2+Z.^2; 

% Decide the working operator
switch parameters.equation
    
    case 'kuprov'
        
        % Isolate rank 2 component of chi
        [~,~,rank2]=mat2sphten(parameters.chi); 
        chi=sphten2mat(0,[0 0 0],rank2);
        
        % Get Kuprov operator in Fourier space
        K=-(1/3)*(chi(1,1)*X.*X+chi(1,2)*X.*Y+chi(1,3)*X.*Z+...
                  chi(2,1)*Y.*X+chi(2,2)*Y.*Y+chi(2,3)*Y.*Z+...
                  chi(3,1)*Z.*X+chi(3,2)*Z.*Y+chi(3,3)*Z.*Z);
       
        % Get Suturina operator in Fourier space
        S=K./L; S(isinf(S)|isnan(S))=0; clear('K','X','Y','Z');
        
    case 'poisson'
        
        % Get the inverse Laplacian in Fourier space
        S=1./L; S(isinf(S)|isnan(S))=0; clear('X','Y','Z');
        
end

% Move operators to the GPU
if parameters.gpu
    S=gpuArray(S); P=gpuArray(P); L=gpuArray(L);
    expt_pcs=gpuArray(parameters.expt_pcs);
else
    expt_pcs=parameters.expt_pcs;
end
              
    % Objective function
    function [err,grad_err,x]=myobj(x)
        
        % Embed the active region
        den=zeros(npoints^3,1,data_type); den(active_point_index)=x;

        % Fourier transform the density
        fft_den=fftn(reshape(den,[npoints npoints npoints]));
        
        % Compute the PCS cube in ppm
        pcs_cube=1e6*real(ifftn(S.*fft_den))/probdens_scaling;
        
        % Get and print cube edge diagnostics
        if ismember('diagnostics',parameters.plot)
            max_abs_edge_val=max([max(max(abs(squeeze(pcs_cube(:,:,1)))))...
                                  max(max(abs(squeeze(pcs_cube(:,:,end)))))...
                                  max(max(abs(squeeze(pcs_cube(:,1,:)))))...
                                  max(max(abs(squeeze(pcs_cube(:,end,:)))))...
                                  max(max(abs(squeeze(pcs_cube(1,:,:)))))...
                                  max(max(abs(squeeze(pcs_cube(end,:,:)))))]);
            disp(['Maximum absolute PCS cube edge value, ppm: ' num2str(max_abs_edge_val)]);
        end
        
        % Project out values at the nuclei
        theo_pcs=P*pcs_cube(:); clear('pcs_cube');
        
        % Compute the least squares error
        err=norm(theo_pcs-expt_pcs,2)^2; err_ls=err;
        
        % Compute contrast functional
        if parameters.sharpen>0
            pivot=max(abs(den))/2; pden=den/pivot; 
            reg_a=contrast_scaling*parameters.sharpen*...
                  sum((pden.^2).*exp(-pden.^2)); 
            err=err+reg_a; clear('pden');
        else
            reg_a=0;
        end
        
        % Compute Tikhonov functional
        if lambda>0
            reg_b=tikhonov_scaling*lambda*...
                  sum(sum(sum(real(ifftn(L.*fft_den)).^2))); 
            err=err+reg_b;
        else
            reg_b=0;
        end
        
        % Compute error gradient
        if nargout > 1
            
            % Fourier transform the density variation
            fftn_dx=fftn(reshape(P.'*(theo_pcs-expt_pcs),[npoints npoints npoints]));
            
            % Get the error gradient
            grad_err=2*1e6*reshape(real(ifftn(S.*fftn_dx)),[npoints^3 1])/probdens_scaling;
            
            % Add contrast gradient (ugly interpolation over the singularity)
            if parameters.sharpen>0
                [pivot,where]=max(abs(den)); pivot=pivot/2; pden=den/pivot;
                grad_cont=2*contrast_scaling*parameters.sharpen*...
                            exp(-pden.^2).*pden.*(1-pden.^2)/pivot;
                grad_cont(where)=(grad_cont(where+1)+grad_cont(where-1))/2; 
                grad_err=grad_err+grad_cont; clear('grad_cont','pden');
            end
            
            % Add Tikhonov gradient
            if lambda>0
                grad_err=grad_err+2*tikhonov_scaling*lambda*...
                         reshape(real(ifftn(L.*L.*fft_den)),[npoints^3 1]);
            end
            
            % Pick out the active region
            grad_err=grad_err(active_point_index);
            
            % Return the gradient to the CPU
            grad_err=gather(grad_err);
                        
        end
        
        % Compute the normalization integral
        normint=(1/probdens_scaling)*(ranges(2)-ranges(1))*...
                                     (ranges(4)-ranges(3))*...
                                     (ranges(6)-ranges(5))*sum(den)/npoints^3;
                 
        % Return the outputs to the CPU
        if parameters.gpu
            err=gather(err); normint=gather(normint); den=gather(den);
            theo_pcs=gather(theo_pcs); expt_pcs=gather(expt_pcs);
            reg_a=gather(reg_a); reg_b=gather(reg_b); err_ls=gather(err_ls);
        end
        
        % Get box corners
        a(1)=parameters.box_cent(1)-parameters.box_size(1)/2; 
        a(2)=parameters.box_cent(1)+parameters.box_size(1)/2;
        a(3)=parameters.box_cent(2)-parameters.box_size(2)/2; 
        a(4)=parameters.box_cent(2)+parameters.box_size(2)/2;
        a(5)=parameters.box_cent(3)-parameters.box_size(3)/2; 
        a(6)=parameters.box_cent(3)+parameters.box_size(3)/2;
        
        % Update the plots
        if ismember('density',parameters.plot)
            set(groot,'CurrentFigure',1);
            volplot(reshape(den,[npoints npoints npoints]),ranges);
        end
        if ismember('box',parameters.plot)
            set(groot,'CurrentFigure',1); hold on;
            plot3([a(1) a(1) a(1) a(1) a(1) a(2) a(2) a(1) a(1) a(2) a(2) a(2) a(2) a(2) a(2) a(1)],...
                  [a(3) a(3) a(4) a(4) a(3) a(3) a(3) a(3) a(4) a(4) a(3) a(3) a(4) a(4) a(4) a(4)],...
                  [a(6) a(5) a(5) a(6) a(6) a(6) a(5) a(5) a(5) a(5) a(5) a(6) a(6) a(5) a(6) a(6)],'r-');
        end
        if ismember('molecule',parameters.plot)
            set(groot,'CurrentFigure',1); hold on; molplot(parameters.xyz_all,conmatrix); kgrid; hold off;
        end
        if ismember('tightzoom',parameters.plot)
            axis([min([parameters.xyz_all(:,1); a(1)]) max([parameters.xyz_all(:,1); a(2)])...
                  min([parameters.xyz_all(:,2); a(3)]) max([parameters.xyz_all(:,2); a(4)])...
                  min([parameters.xyz_all(:,3); a(5)]) max([parameters.xyz_all(:,3); a(6)])]);
        end
        if ismember('diagnostics',parameters.plot)
            set(groot,'CurrentFigure',2); clf reset; hold on; plot(expt_pcs,theo_pcs,'ro');
            plot([min(expt_pcs) max(expt_pcs)],[min(expt_pcs) max(expt_pcs)],'b-');
            annotation('textbox',[0.15 0.7 0.2 0.2],...
                       'String',{['Density integral: '  num2str(normint)],...
                                 ['Least sq. penalty: ' num2str(err_ls)],...
                                 ['Contrast penalty: '  num2str(reg_a)],...
                                 ['Tikhonov penalty: '  num2str(reg_b)],...
                                 ['Total penalty: '     num2str(err)]},'FitBoxToText','on');
            kxlabel('Experimental PCS, ppm'); kylabel('Predicted PCS, ppm'); kgrid; box on;
        end
        drawnow();
              
    end

    % Hessian-times-vector function
    function y=myhess(x,v)
        
        % Preallocate the answer
        y=zeros(npoints^3,size(v,2),data_type);
        
        % Loop over the vectors in the answer
        for n=1:size(v,2)
            
            % Preallocate arrays
            den_v=zeros(npoints^3,1,data_type); 
            den_x=zeros(npoints^3,1,data_type); 
            
            % Embed the active region
            den_v(active_point_index)=v(:,n);
            den_x(active_point_index)=x;
            
            % Run the first FFT
            fftn_probden=fftn(reshape(den_v,[npoints npoints npoints]));
            theo_pcs=1e6*real(ifftn(S.*fftn_probden))/probdens_scaling;
        
            % Filter through nuclear positions
            probden=reshape(P.'*(P*theo_pcs(:)),[npoints npoints npoints]);
            
            % Run the second FFT
            y(:,n)=2*1e6*reshape(real(ifftn(S.*fftn(probden))),[npoints^3 1])/probdens_scaling;
            
            % Add contrast Hessian (ugly interpolation over the singularity)
            if parameters.sharpen>0
                [pivot,where]=max(abs(den_x)); pivot=pivot/2; pden_x=den_x/pivot;
                hess_cont=2*contrast_scaling*parameters.sharpen*(1/pivot^2)*...
                            exp(-pden_x.^2).*(1-5*pden_x.^2+2*pden_x.^4).*den_v;
                hess_cont(where)=(hess_cont(where+1)+hess_cont(where-1))/2;
                y(:,n)=y(:,n)+hess_cont;
            end
            
            % Add Tikhonov Hessian
            if lambda>0
                y(:,n)=y(:,n)+2*tikhonov_scaling*lambda*...
                              reshape(real(ifftn(L.*L.*fftn_probden)),[npoints^3 1]);
            end
            
        end
        
        % Pick out the active region
        y=y(active_point_index,:);
        
        % Return to the CPU
        if parameters.gpu, y=gather(y); end
        
    end

% Get the figures going
figure(2); figure(1);

% Set the optimiser to full Newton with implicit Hessian
options=optimoptions('fmincon','Algorithm','trust-region-reflective','Display','iter','GradObj','on',...
                     'Hessian','user-supplied','MaxIterations',Inf,'FunctionTolerance',1e-12,...
                     'OptimalityTolerance',1e-12,'StepTolerance',1e-12,'MaxFunctionEvaluations',Inf,...
                     'HessMult',@myhess,'SubproblemAlgorithm','factorization');

% Decide the boundaries
switch parameters.equation
    
    case 'kuprov'
        
        % Bounded from below
        lower_bound=zeros(size(guess)); upper_bound=[];
        
    case 'poisson'
        
        % Unbounded
        lower_bound=[]; upper_bound=[];
        
end
                 
% Run the error minimization
answer=fmincon(@myobj,guess,[],[],[],[],lower_bound,upper_bound,[],options);

% Embed the active region
source_cube=zeros([npoints npoints npoints]);
source_cube(active_point_index)=answer/probdens_scaling;

% Back-calculate pseudocontact shifts
pcs_cube=1e6*real(ifftn(S.*fftn(source_cube))); pred_pcs=gather(P*pcs_cube(:));

end

% Consistency enforcement
function grumble(parameters,npoints,lambda)
if verLessThan('matlab','9.0')
    error('3D PCS reconstruction module requires 64-bit Matlab R2016a or later.');
end
if verLessThan('parallel','6.8')
    error('3D PCS reconstruction module requires Matlab Parallel Computing Toolbox version 6.8 or later.');
end
if verLessThan('optim','7.4')
    error('3D PCS reconstruction module requires Matlab Optimisation Toolbox version 7.4 or later.');
end
if ~isfield(parameters,'equation')
    error('working equation must be specified in parameters.equation field.');
end
if ~ismember(parameters.equation,{'kuprov','poisson'})
    error('parameters.equation may be set to ''kuprov'' and ''poisson''.');
end
if ~isfield(parameters,'xyz')
    error('coordinates of PCS active nuclei must be provided in parameters.xyz field.');
elseif (~isnumeric(parameters.xyz))||(~isreal(parameters.xyz))||(size(parameters.xyz,2)~=3)
    error('parameters.xyz must be an Nx3 array of atomic coordinates.');
end
if ~isfield(parameters,'expt_pcs')
    error('experimental PCS values must be provided in parameters.expt_pcs field.');
elseif (~isnumeric(parameters.expt_pcs))||(size(parameters.expt_pcs,2)~=1)||(~isreal(parameters.expt_pcs))
    error('parameters.expt_pcs must be a real column vector.');
end
if ~isfield(parameters,'xyz_all')
    error('coordinates of all atoms in the molecule must be provided in parameters.xyz_all field.');
elseif (~isnumeric(parameters.xyz_all))||(~isreal(parameters.xyz_all))||(size(parameters.xyz_all,2)~=3)
    error('parameters.xyz_all must be an Nx3 array of atomic coordinates.');
end
if size(parameters.xyz,1)~=size(parameters.expt_pcs,1)
    error('the number of rows in parameters.xyz and parameters.expt_pcs must be the same.');
end
if ~isfield(parameters,'plot')
    error('visualisation options must be supplied in parameters.plot field.');
elseif (~iscell(parameters.plot))||any(~cellfun(@ischar,parameters.plot))
    error('all elements of parameters.plot cell array must be character strings.');
elseif ~isempty(setdiff(parameters.plot,{'diagnostics','density','molecule','tightzoom','box'}))
    error('valid elements of parameters.plot are ''diagnostics'', ''density'', ''molecule'', ''tightzoom'', and ''box''.');
end
if ~isfield(parameters,'box_cent')
    error('source box centre coordinates must be supplied in parameters.box_cent field.');
elseif (~isnumeric(parameters.box_cent))||(~isreal(parameters.box_cent))||(numel(parameters.box_cent)~=3)
    error('parameters.box_cent must be a real vector with three elements.');
end
if ~isfield(parameters,'box_size')
    error('source box sizes must be supplied in parameters.box_size field.');
elseif (~isnumeric(parameters.box_size))||(~isreal(parameters.box_size))||(numel(parameters.box_size)~=3)
    error('parameters.box_size must be a real vector with three elements.');
end
if ~isfield(parameters,'margins')
    error('solution box margin widths must be supplied in parameters.margins field.');
elseif (~isnumeric(parameters.margins))||(~isreal(parameters.margins))||(numel(parameters.margins)~=6)
    error('parameters.margins must be a real vector with six elements.');
end
if ~isfield(parameters,'confine')
    error('source confinement radii must be supplied in parameters.confine field.');
elseif (~isnumeric(parameters.confine))||(~isreal(parameters.confine))||...
       (numel(parameters.confine)~=2)||any(parameters.confine<0)
    error('parameters.confine must be a real vector with two non-negative elements.');
end
if ~isfield(parameters,'sharpen')
    error('contrast parameter must be supplied in parameters.sharpen field.');
elseif (~isnumeric(parameters.sharpen))||(~isreal(parameters.sharpen))||(numel(parameters.sharpen)~=1)
    error('parameters.sharpen must be a real scalar.');
end
if ~isfield(parameters,'chi')
    error('effective susceptibility tensor must be supplied in parameters.chi field.');
elseif (~isnumeric(parameters.chi))||(~isreal(parameters.chi))||(any(size(parameters.chi)~=3))
    error('parameters.chi must be a real symmetric 3x3 matrix.');
end
if (~isnumeric(npoints))||(~isscalar(npoints))||(mod(npoints,1)~=0)||(npoints<=10)
    error('npoints must be a positive integer greater than 10.');
end
if (~isnumeric(lambda))||(~isreal(lambda))
    error('lambda must be a real vector.');
end
end

% People call me a perfectionist, but I'm not. I'm a
% rightist. I do something until it's right, and then
% I move on to the next thing.
%
% James Cameron

