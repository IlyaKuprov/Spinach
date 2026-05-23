% Adaptively recursed Voitlander integrator. Computes an approximation
% of an integral of field-swept EPR transition over a spherical triang-
% le. Syntax:
%
%             spec=voitlander(spin_system,parameters,...
%                             triangle,Ic,Iz,Qc,Qz,Hmw)
%
% Parameters:
%
%     triangle(1:3).xyz - Cartesian coordinates of the corners of the
%                         spherical triangle, unit column vectors
%
%     triangle(1:3).tf  - transition fields at the corners of the
%                         spherical triangle, real column vectors, one
%                         element per transition
%
%     triangle(1:3).tm  - transition moments at the corners of the sphe-
%                         rical triangle, positive column vectors, one
%                         element per transition
%
%     triangle(1:3).tw  - transition widths at the corners of the sphe-
%                         rical triangle, positive column vectors, one
%                         element per transition
%
%     triangle(1:3).pd  - energy level population differences at the
%                         triangle corners, real column vectors, one
%                         element per transition
%
%     triangle(1:3).ti  - transition identity arrays at the triangle
%                         corners, one row per transition
%
%     triangle(1:3).tj  - scaled field-sweep Jacobians at the triangle
%                         corners, real column vectors, one element per
%                         transition
%
%     Ic          - isotropic part of the coupling Hamiltonian,
%                   a Hermitian matrix (set retention to 'couplings'
%                   in assume.m and then call hamiltonian.m)
%
%     Qc          - irreducible components of the anisotropic part
%                   of the coupling Hamiltonian a cell array re-
%                   turned by hamiltonian.m
%
%     Iz          - isotropic part of the Zeeman Hamiltonian, a
%                   Hermitian matrix (set retention to 'zeeman'
%                   in assume.m and then call hamiltonian.m) nor-
%                   malised to 1 Tesla
%
%     Qz          - irreducible components of the anisotropic part
%                   of the Zeeman Hamiltonian a cell array retur-
%                   ned by hamiltonian.m, normalised to 1 Tesla
%
%     Hmw         - perturbation operator, a Hermitian matrix
%
%     parameters.b_axis  - a vector of magnetic field values, Tesla
%
%     parameters.int_tol - integration accuracy tolerance
%
% Outputs:
%
%     spec        - ESR spectrum integral over the triangle, array
%                   of the same dimension as parameters.b_axis
%               
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=voitlander.m>

function spec=voitlander(spin_system,parameters,triangle,Ic,Iz,Qc,Qz,Hmw)

% Check consistency
grumble(spin_system,parameters,triangle,Ic,Iz,Qc,Qz,Hmw);

% Get traingle subdivision midpoints
[r12,r23,r31]=sphtrsubd(triangle(1).xyz,...
                        triangle(2).xyz,...
                        triangle(3).xyz);

% First subdivision midpoint eigenset
[phi,elev]=cart2sph(r12(1),r12(2),r12(3));
parameters12=parameters; 
parameters12.orientation=[0, pi/2-elev, phi];
Hz=Iz+orientation(Qz,parameters12.orientation); 
Hc=Ic+orientation(Qc,parameters12.orientation); 
eigenset12=eigenfields(spin_system,parameters12,Hz,Hc,Hmw);
eigenset12.xyz=[r12(1); r12(2); r12(3)];

% Second subdivision midpoint eigenset
[phi,elev]=cart2sph(r23(1),r23(2),r23(3));
parameters23=parameters; 
parameters23.orientation=[0, pi/2-elev, phi];
Hz=Iz+orientation(Qz,parameters23.orientation); 
Hc=Ic+orientation(Qc,parameters23.orientation); 
eigenset23=eigenfields(spin_system,parameters23,Hz,Hc,Hmw);
eigenset23.xyz=[r23(1); r23(2); r23(3)];

% Third subdivision midpoint eigenset
[phi,elev]=cart2sph(r31(1),r31(2),r31(3));
parameters31=parameters; 
parameters31.orientation=[0, pi/2-elev, phi];
Hz=Iz+orientation(Qz,parameters31.orientation); 
Hc=Ic+orientation(Qc,parameters31.orientation); 
eigenset31=eigenfields(spin_system,parameters31,Hz,Hc,Hmw);
eigenset31.xyz=[r31(1); r31(2); r31(3)];

% Eigensets for the vertices (three each) of the four subdivision traingles
triangle_a(1)=triangle(1); triangle_a(2)=eigenset12;  triangle_a(3)=eigenset31;
triangle_b(1)=eigenset12;  triangle_b(2)=triangle(2); triangle_b(3)=eigenset23;
triangle_c(1)=eigenset31;  triangle_c(2)=eigenset23;  triangle_c(3)=triangle(3);
triangle_d(1)=eigenset12;  triangle_d(2)=eigenset23;  triangle_d(3)=eigenset31;

% Compute the subdivided integral
spec_sub=trint(parameters,triangle_a)+...
         trint(parameters,triangle_b)+...
         trint(parameters,triangle_c)+...
         trint(parameters,triangle_d);
     
% Compute the direct integral
spec_dir=trint(parameters,triangle);

% If the accuracy is insufficient, recurse
if norm(spec_dir-spec_sub,2)>parameters.int_tol

    % Compute the four triangles of the subdivision via asynchronous recursion
    spec_a=parfeval(@voitlander,1,spin_system,parameters,triangle_a,Ic,Iz,Qc,Qz,Hmw);
    spec_b=parfeval(@voitlander,1,spin_system,parameters,triangle_b,Ic,Iz,Qc,Qz,Hmw);
    spec_c=parfeval(@voitlander,1,spin_system,parameters,triangle_c,Ic,Iz,Qc,Qz,Hmw);
    spec_d=parfeval(@voitlander,1,spin_system,parameters,triangle_d,Ic,Iz,Qc,Qz,Hmw);

    % Retrieve the work
    spec=fetchOutputs(spec_a)+fetchOutputs(spec_b)+...
         fetchOutputs(spec_c)+fetchOutputs(spec_d);

else
    
    % Original triangle
    spec=spec_sub;

end

end

% Single-triangle Voitlander integrator
function spec=trint(parameters,triangle)
                 
% Preallocate the spectrum
spec=zeros(size(parameters.b_axis),'like',1i);

% Skip missing corners
if isempty(triangle(1).ti)||...
   isempty(triangle(2).ti)||...
   isempty(triangle(3).ti)
    return;
end

% Initialise transition index lists
idx1=[]; idx2=[]; idx3=[];
used1=false(numel(triangle(1).tf),1);
used2=false(numel(triangle(2).tf),1);
used3=false(numel(triangle(3).tf),1);

% Get spherical edge lengths
edge12=max(acos(max(-1,min(1,dot(triangle(1).xyz(:),...
                                      triangle(2).xyz(:))))),eps);
edge23=max(acos(max(-1,min(1,dot(triangle(2).xyz(:),...
                                      triangle(3).xyz(:))))),eps);
edge31=max(acos(max(-1,min(1,dot(triangle(3).xyz(:),...
                                      triangle(1).xyz(:))))),eps);

% Process level pairs
if (size(triangle(1).ti,2)>=2)&&...
   (size(triangle(2).ti,2)>=2)&&...
   (size(triangle(3).ti,2)>=2)

    % Find level pairs with roots at triangle vertices
    pair_list=unique([triangle(1).ti(:,1:2); ...
                      triangle(2).ti(:,1:2); ...
                      triangle(3).ti(:,1:2)],'rows','stable');

    % Match same-pair roots by
    % nearest field continuation
    for p=1:size(pair_list,1)

        % Get candidate roots for this level pair
        cand1=find((~used1)&...
                   (triangle(1).ti(:,1)==pair_list(p,1))&...
                   (triangle(1).ti(:,2)==pair_list(p,2)));
        cand2=find((~used2)&...
                   (triangle(2).ti(:,1)==pair_list(p,1))&...
                   (triangle(2).ti(:,2)==pair_list(p,2)));
        cand3=find((~used3)&...
                   (triangle(3).ti(:,1)==pair_list(p,1))&...
                   (triangle(3).ti(:,2)==pair_list(p,2)));
        if isempty(cand1)||isempty(cand2)||isempty(cand3), continue; end

        % Match roots by global field-continuation cost
        [root1,root2,root3]=rootmatch(triangle(1).tf(cand1),...
                                      triangle(2).tf(cand2),...
                                      triangle(3).tf(cand3),...
                                      edge12,edge23,edge31);

        % Index updates
        idx1=[idx1 cand1(root1).']; %#ok<AGROW>
        idx2=[idx2 cand2(root2).']; %#ok<AGROW>
        idx3=[idx3 cand3(root3).']; %#ok<AGROW>
        used1(cand1(root1))=true;
        used2(cand2(root2))=true;
        used3(cand3(root3))=true;

    end

end

% Match identity-poor leftovers by field continuation
if (size(triangle(1).ti,2)<2)||...
   (size(triangle(2).ti,2)<2)||...
   (size(triangle(3).ti,2)<2)

    % Get unused roots
    cand1=find(~used1);
    cand2=find(~used2);
    cand3=find(~used3);

    % Match roots by global field-continuation cost
    [root1,root2,root3]=rootmatch(triangle(1).tf(cand1),...
                                  triangle(2).tf(cand2),...
                                  triangle(3).tf(cand3),...
                                  edge12,edge23,edge31);

    % Index updates
    idx1=[idx1 cand1(root1).'];
    idx2=[idx2 cand2(root2).'];
    idx3=[idx3 cand3(root3).'];

end

% Transition count
ntrans=numel(idx1);

% Get triangle area
S=sphtarea(triangle(1).xyz,...
           triangle(2).xyz,...
           triangle(3).xyz);

% Build spectrum
for n=1:ntrans

    % Get transition indices
    idx_a=idx1(n); idx_b=idx2(n); idx_c=idx3(n);

    % Get signal information
    line_width=mean([triangle(1).tw(idx_a) ...
                     triangle(2).tw(idx_b) ...
                     triangle(3).tw(idx_c)]);

    % Average finite vertex-wise intensity weights
    tran_prod=[triangle(1).tm(idx_a)*triangle(1).tj(idx_a)*triangle(1).pd(idx_a) ...
               triangle(2).tm(idx_b)*triangle(2).tj(idx_b)*triangle(2).pd(idx_b) ...
               triangle(3).tm(idx_c)*triangle(3).tj(idx_c)*triangle(3).pd(idx_c)];
    tran_prod=tran_prod(isfinite(tran_prod));

    % Drop instabilities
    if isempty(tran_prod)
        continue;
    else
        tran_amp=mean(tran_prod);
    end

    % Convolutions of Lorentzians with triangles
    spec=spec+S*lorentzcon([triangle(1).tf(idx_a) ...
                            triangle(2).tf(idx_b)...
                            triangle(3).tf(idx_c)],tran_amp,...
                            line_width,parameters.b_axis);
    
end

end

% Consistency enforcement
function grumble(spin_system,parameters,triangle,Ic,Iz,Qc,Qz,Hmw)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||...
   (~isfield(spin_system.bas,'formalism'))
    error('spin_system must be a Spinach data structure with basis information.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isfield(parameters,'b_axis')
    error('field sweep axis must be supplied in parameters.b_axis field.');
end
if (~isnumeric(parameters.b_axis))||(~isreal(parameters.b_axis))||...
   (~isvector(parameters.b_axis))||any(~isfinite(parameters.b_axis(:)))
    error('parameters.b_axis must be a finite real vector.');
end
if ~isfield(parameters,'int_tol')
    error('integration tolerance must be supplied in parameters.int_tol field.');
end
if (~isnumeric(parameters.int_tol))||(~isreal(parameters.int_tol))||...
   (~isscalar(parameters.int_tol))||(parameters.int_tol<=0)
    error('parameters.int_tol must be a positive real scalar.');
end
if ~isfield(parameters,'mw_freq')
    error('resonance frequency must be supplied in parameters.mw_freq field.');
end
if (~isnumeric(parameters.mw_freq))||(~isreal(parameters.mw_freq))||...
   (~isscalar(parameters.mw_freq))
    error('parameters.mw_freq must be a real scalar.');
end
if ~isfield(parameters,'window')
    error('field window must be supplied in parameters.window field.');
end
if (~isnumeric(parameters.window))||(~isreal(parameters.window))||...
   (numel(parameters.window)~=2)
    error('parameters.window must have two real elements.');
end
if ~isfield(parameters,'pp_tol')
    error('peak position tolerance must be supplied in parameters.pp_tol field.');
end
if (~isnumeric(parameters.pp_tol))||(~isreal(parameters.pp_tol))||...
   (~isscalar(parameters.pp_tol))
    error('parameters.pp_tol must be a real scalar.');
end
if ~isfield(parameters,'tm_tol')
    error('transition moment tolerance must be supplied in parameters.tm_tol field.');
end
if (~isnumeric(parameters.tm_tol))||(~isreal(parameters.tm_tol))||...
   (~isscalar(parameters.tm_tol))
    error('parameters.tm_tol must be a real scalar.');
end
if ~isfield(parameters,'fwhm')
    error('transition FWHM must be supplied in parameters.fwhm field.');
end
if (~isnumeric(parameters.fwhm))||(~isreal(parameters.fwhm))||...
   (~isscalar(parameters.fwhm))||(~isfinite(parameters.fwhm))||...
   (parameters.fwhm<=0)
    error('parameters.fwhm must be a finite positive real scalar.');
end
if strcmp(spin_system.bas.formalism,'zeeman-hilb')
    if ~isfield(parameters,'rspt_order')
        error('perturbation theory order must be supplied in parameters.rspt_order field.');
    end
    if (~isnumeric(parameters.rspt_order))||(~isreal(parameters.rspt_order))||...
       (~isscalar(parameters.rspt_order))||((mod(parameters.rspt_order,1)~=0)&&...
       (~isinf(parameters.rspt_order)))||(parameters.rspt_order<0)
        error('parameters.rspt_order must be a non-negative integer or Inf.');
    end
end
if (~isstruct(triangle))||(numel(triangle)~=3)
    error('triangle must be a three-element eigenset structure array.');
end
if ~all(isfield(triangle,{'xyz','tf','tm','tw','pd','ti','tj'}))
    error('triangle vertices must contain xyz, tf, tm, tw, pd, ti, and tj fields.');
end
for n=1:3
    if (~isnumeric(triangle(n).xyz))||(~isreal(triangle(n).xyz))||...
       (numel(triangle(n).xyz)~=3)||any(~isfinite(triangle(n).xyz(:)))
        error('triangle vertex xyz fields must be three-element finite real vectors.');
    end
    if abs(norm(triangle(n).xyz,2)-1)>sqrt(eps)
        error('triangle vertex xyz fields must be unit vectors.');
    end
    if (~isnumeric(triangle(n).tf))||(~isreal(triangle(n).tf))||...
       (size(triangle(n).tf,2)~=1)||any(~isfinite(triangle(n).tf(:)))
        error('triangle vertex tf fields must be finite real column vectors.');
    end
    if (~isnumeric(triangle(n).tm))||(~isreal(triangle(n).tm))||...
       (size(triangle(n).tm,2)~=1)||any(~isfinite(triangle(n).tm(:)))||...
       any(triangle(n).tm(:)<0)
        error('triangle vertex tm fields must be finite non-negative real column vectors.');
    end
    if (~isnumeric(triangle(n).tw))||(~isreal(triangle(n).tw))||...
       (size(triangle(n).tw,2)~=1)||any(~isfinite(triangle(n).tw(:)))||...
       any(triangle(n).tw(:)<=0)
        error('triangle vertex tw fields must be finite positive real column vectors.');
    end
    if (~isnumeric(triangle(n).pd))||(~isreal(triangle(n).pd))||...
       (size(triangle(n).pd,2)~=1)||any(~isfinite(triangle(n).pd(:)))
        error('triangle vertex pd fields must be finite real column vectors.');
    end
    if (~isnumeric(triangle(n).tj))||(~isreal(triangle(n).tj))||...
       (size(triangle(n).tj,2)~=1)||any(isnan(triangle(n).tj(:)))
        error('triangle vertex tj fields must be real column vectors without NaN elements.');
    end
    if (~isnumeric(triangle(n).ti))||(~isreal(triangle(n).ti))||...
       any(~isfinite(triangle(n).ti(:)))||any(mod(triangle(n).ti(:),1)~=0)||...
       any(triangle(n).ti(:)<1)
        error('triangle vertex ti fields must contain finite positive real integers.');
    end
    if (numel(triangle(n).tf)~=numel(triangle(n).tm))||...
       (numel(triangle(n).tf)~=numel(triangle(n).tw))||...
       (numel(triangle(n).tf)~=numel(triangle(n).pd))||...
       (numel(triangle(n).tf)~=numel(triangle(n).tj))||...
       (numel(triangle(n).tf)~=size(triangle(n).ti,1))
        error('triangle vertex transition fields must have matching numbers of rows.');
    end
end
if (~isnumeric(Ic))||(~isnumeric(Iz))||(~ismatrix(Ic))||...
   (~ismatrix(Iz))||(size(Ic,1)~=size(Ic,2))||...
   (size(Iz,1)~=size(Iz,2))||(~all(size(Ic)==size(Iz)))
    error('Ic and Iz must be square matrices of the same size.');
end
if (~ishermitian(Ic))||(~ishermitian(Iz))
    error('Ic and Iz must be Hermitian.');
end
if (~iscell(Qc))||(~iscell(Qz))
    error('Qc and Qz must be cell arrays.');
end
if (~isnumeric(Hmw))||(~ismatrix(Hmw))||...
   ((~isequal(size(Hmw),size(Ic)))&&(~isequal(size(Hmw),[size(Ic,1) 1])))
    error('Hmw must be a matrix matching Ic and Iz or a column vector of matching dimension.');
end

end

% "How many divisions?"
%
% Joseph Stalin, when told that
% Vatican was a powerful entity
