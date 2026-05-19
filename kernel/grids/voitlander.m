% Adaptively recursed Voitlander integrator. Computes an approxi-
% mation of an integral of field-swept EPR transition over a sph-
% erical triangle. Syntax:
%
% spec=voitlander(spin_system,parameters,tri,Ic,Iz,Qc,Qz,Hmw,b_axis)
%
% Parameters:
%
%     tri.vert(n).xyz - Cartesian coordinates of the corners of the
%                       spherical triangle, unit column vectors
%
%     tri.vert(n).tf  - transition fields at the corners of the
%                       spherical triangle, real column vectors, one
%                       element per transition
%
%     tri.vert(n).tm  - transition moments at the corners of the sphe-
%                       rical triangle, positive column vectors, one
%                       element per transition
%
%     tri.vert(n).tw  - transition widths at the corners of the sphe-
%                       rical triangle, positive column vectors, one
%                       element per transition
%
%     tri.vert(n).pd  - energy level population differences at the
%                       triangle corners, real column vectors, one
%                       element per transition
%
%     tri.vert(n).ti  - transition identity arrays at the triangle
%                       corners, one row per transition
%
%     tri.vert(n).tj  - scaled field-sweep Jacobians at the triangle
%                       corners, real column vectors, one element per
%                       transition
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
%     b_axis      - a vector of magnetic field values, Tesla
%
%     parameters.int_tol - recursive integration accuracy
%                          tolerance
%
% Outputs:
%
%     spec        - ESR spectrum integral over the triangle, array
%                   of the same dimension as b_axis
%               
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=voitlander.m>

function spec=voitlander(spin_system,parameters,tri,Ic,Iz,Qc,Qz,Hmw,b_axis)

% Check consistency
grumble(parameters,tri,Ic,Iz,Qc,Qz,b_axis);

% Extract triangle vertices
vert1=tri.vert(1); vert2=tri.vert(2); vert3=tri.vert(3);

% Subdivide the triangle
[r12,r23,r31]=sphtrsubd(vert1.xyz,vert2.xyz,vert3.xyz);

% Characterise new vertices
vert12=struct(); vert12.xyz=r12;
[phi,elev]=cart2sph(r12(1),r12(2),r12(3));
theta=pi/2-elev; loc_params=parameters; loc_params.orientation=[0 theta phi];
[vert12.tf,vert12.tm,vert12.tw,vert12.pd,vert12.ti,vert12.tj]=...
    eigenfields(spin_system,loc_params,Iz,Qz,Ic,Qc,Hmw);

vert23=struct(); vert23.xyz=r23;
[phi,elev]=cart2sph(r23(1),r23(2),r23(3)); 
theta=pi/2-elev; loc_params=parameters; loc_params.orientation=[0 theta phi];
[vert23.tf,vert23.tm,vert23.tw,vert23.pd,vert23.ti,vert23.tj]=...
    eigenfields(spin_system,loc_params,Iz,Qz,Ic,Qc,Hmw);

vert31=struct(); vert31.xyz=r31;
[phi,elev]=cart2sph(r31(1),r31(2),r31(3));
theta=pi/2-elev; loc_params=parameters; loc_params.orientation=[0 theta phi];
[vert31.tf,vert31.tm,vert31.tw,vert31.pd,vert31.ti,vert31.tj]=...
    eigenfields(spin_system,loc_params,Iz,Qz,Ic,Qc,Hmw);

% Assemble the subdivided triangles
tri_a.vert=[vert1 vert12 vert31];
tri_b.vert=[vert12 vert2 vert23];
tri_c.vert=[vert31 vert23 vert3];
tri_d.vert=[vert12 vert23 vert31];

% Compute the subdivided integral
spec_sub=trint(tri_a,b_axis)+trint(tri_b,b_axis)+...
         trint(tri_c,b_axis)+trint(tri_d,b_axis);
     
% Compute the direct integral
spec_dir=trint(tri,b_axis);

% If the accuracy is insufficient, make a recursive call
if norm(spec_dir-spec_sub,2)>parameters.int_tol

    % Four triangles of the subdivision
    spec=voitlander(spin_system,parameters,tri_a,Ic,Iz,Qc,Qz,Hmw,b_axis)+...
         voitlander(spin_system,parameters,tri_b,Ic,Iz,Qc,Qz,Hmw,b_axis)+...
         voitlander(spin_system,parameters,tri_c,Ic,Iz,Qc,Qz,Hmw,b_axis)+...
         voitlander(spin_system,parameters,tri_d,Ic,Iz,Qc,Qz,Hmw,b_axis);
else
    
    % Original triangle
    spec=spec_sub;

end

end

% Single triangle integrator
function spec=trint(tri,b_axis)
                 
% Preallocate the spectrum
spec=zeros(size(b_axis),'like',1i);

% Extract triangle vertices
vert1=tri.vert(1); vert2=tri.vert(2); vert3=tri.vert(3);

% Find transitions present at all three vertices
if isempty(vert1.ti)||isempty(vert2.ti)||isempty(vert3.ti), return; end
if size(vert1.ti,2)>=3

    % Find level pairs with roots at triangle vertices
    pair_list=unique([vert1.ti(:,1:2); vert2.ti(:,1:2);...
                      vert3.ti(:,1:2)],'rows','stable');
    idx1=[]; idx2=[]; idx3=[];

    % Match same-pair roots by nearest field continuation
    for p=1:size(pair_list,1)

        % Get candidate roots for this level pair
        cand1=find((vert1.ti(:,1)==pair_list(p,1))&(vert1.ti(:,2)==pair_list(p,2)));
        cand2=find((vert2.ti(:,1)==pair_list(p,1))&(vert2.ti(:,2)==pair_list(p,2)));
        cand3=find((vert3.ti(:,1)==pair_list(p,1))&(vert3.ti(:,2)==pair_list(p,2)));
        if isempty(cand1)||isempty(cand2)||isempty(cand3), continue; end

        % Initialise local candidate masks
        used1=false(size(cand1)); used2=false(size(cand2)); used3=false(size(cand3));

        % Greedily match roots with the smallest field spread
        while any(~used1)&&any(~used2)&&any(~used3)
            best_span=Inf; best_pick=[0 0 0];
            list1=find(~used1).'; list2=find(~used2).'; list3=find(~used3).';
            for a=list1
                for b=list2
                    for c=list3
                        field_vals=[vert1.tf(cand1(a)) vert2.tf(cand2(b))...
                                    vert3.tf(cand3(c))];
                        field_span=max(field_vals)-min(field_vals);
                        if field_span<best_span
                            best_span=field_span; best_pick=[a b c];
                        end
                    end
                end
            end
            if ~isfinite(best_span), break; end
            idx1(end+1)=cand1(best_pick(1)); %#ok<AGROW>
            idx2(end+1)=cand2(best_pick(2)); %#ok<AGROW>
            idx3(end+1)=cand3(best_pick(3)); %#ok<AGROW>
            used1(best_pick(1))=true;
            used2(best_pick(2))=true;
            used3(best_pick(3))=true;
        end
    end
    ntrans=numel(idx1);
else
    [ti12,idx1,idx2]=intersect(vert1.ti,vert2.ti,'rows','stable');
    [~,idx12,idx3]=intersect(ti12,vert3.ti,'rows','stable');
    idx1=idx1(idx12); idx2=idx2(idx12);
    ntrans=numel(idx1);
end

% Area of spherical triangle
S=sphtarea(vert1.xyz,vert2.xyz,vert3.xyz);

% Build spectrum
for n=1:ntrans

    % Get transition indices
    idx_a=idx1(n); idx_b=idx2(n); idx_c=idx3(n);

    % Get signal information
    min_freq=min([vert1.tf(idx_a) vert2.tf(idx_b) vert3.tf(idx_c)]);
    max_freq=max([vert1.tf(idx_a) vert2.tf(idx_b) vert3.tf(idx_c)]);
    tran_mom=mean([vert1.tm(idx_a) vert2.tm(idx_b) vert3.tm(idx_c)]);
    pop_diff=mean([vert1.pd(idx_a) vert2.pd(idx_b) vert3.pd(idx_c)]);
    line_width=mean([vert1.tw(idx_a) vert2.tw(idx_b) vert3.tw(idx_c)]);
    jac_weight=mean([vert1.tj(idx_a) vert2.tj(idx_b) vert3.tj(idx_c)]);

    % Find the relevant part of the axis
    b_mask=(b_axis>min_freq-3*line_width)&(b_axis<max_freq+3*line_width);

    % Convolutions of Lorentzians with triangles
    spec(b_mask)=spec(b_mask)+...
                 S*pop_diff*lorentzcon([vert1.tf(idx_a) vert2.tf(idx_b)...
                                         vert3.tf(idx_c)],tran_mom*jac_weight,...
                                        line_width,b_axis(b_mask));
                  
end

end

% Consistency enforcement
function grumble(parameters,tri,Ic,Iz,Qc,Qz,b_axis)
if (~isstruct(tri))||(~isfield(tri,'vert'))||(numel(tri.vert)~=3)
    error('tri.vert must be a three-element vertex structure array.');
end
vert_fields={'xyz','tf','tm','tw','pd','ti','tj'};
for n=1:numel(vert_fields)
    if ~isfield(tri.vert,vert_fields{n})
        error('tri.vert entries must contain xyz, tf, tm, tw, pd, ti, and tj fields.');
    end
end
for n=1:3
    if (~isnumeric(tri.vert(n).xyz))||(~isreal(tri.vert(n).xyz))||...
       (~iscolumn(tri.vert(n).xyz))||(numel(tri.vert(n).xyz)~=3)||...
       (abs(norm(tri.vert(n).xyz,2)-1)>1e-6)
        error('tri.vert(n).xyz must be a real unit column 3-vector.');
    end
    if (~isnumeric(tri.vert(n).tf))||(~isreal(tri.vert(n).tf))||...
       (~iscolumn(tri.vert(n).tf))
        error('tri.vert(n).tf must be a real column vector.');
    end
    if (~isnumeric(tri.vert(n).tm))||(~isreal(tri.vert(n).tm))||...
       (~iscolumn(tri.vert(n).tm))
        error('tri.vert(n).tm must be a real column vector.');
    end
    if (~isnumeric(tri.vert(n).tw))||(~isreal(tri.vert(n).tw))||...
       (~iscolumn(tri.vert(n).tw))
        error('tri.vert(n).tw must be a real column vector.');
    end
    if (~isnumeric(tri.vert(n).pd))||(~isreal(tri.vert(n).pd))||...
       (~iscolumn(tri.vert(n).pd))
        error('tri.vert(n).pd must be a real column vector.');
    end
    if (~isnumeric(tri.vert(n).ti))||(~isreal(tri.vert(n).ti))||...
       (size(tri.vert(n).ti,1)~=numel(tri.vert(n).tf))
        error('tri.vert(n).ti must be a real array with one row per transition.');
    end
    if (~isnumeric(tri.vert(n).tj))||(~isreal(tri.vert(n).tj))||...
       (~iscolumn(tri.vert(n).tj))
        error('tri.vert(n).tj must be a real column vector.');
    end
    if (numel(tri.vert(n).tm)~=numel(tri.vert(n).tf))||...
       (numel(tri.vert(n).tw)~=numel(tri.vert(n).tf))||...
       (numel(tri.vert(n).pd)~=numel(tri.vert(n).tf))||...
       (numel(tri.vert(n).tj)~=numel(tri.vert(n).tf))
        error('tri.vert(n).tf, tm, tw, pd, ti, and tj sizes must agree.');
    end
end
if (size(tri.vert(1).ti,2)~=size(tri.vert(2).ti,2))||...
   (size(tri.vert(1).ti,2)~=size(tri.vert(3).ti,2))
    error('tri.vert(n).ti arrays must have the same number of columns.');
end
if ~isfield(parameters,'int_tol')
    error('integration accuracy tolerance must be supplied in parameters.int_tol field.');
end
if (~isnumeric(parameters.int_tol))||(~isreal(parameters.int_tol))||...
   (~isscalar(parameters.int_tol))||(parameters.int_tol<=0)
    error('parameters.int_tol must be a positive real scalar.');
end
if (~isnumeric(Ic))||(size(Ic,1)~=size(Ic,2))||...
   (~isnumeric(Iz))||(size(Iz,1)~=size(Iz,2))
    error('Ic and Iz must be square Hermitian matrices.');
end
if (~iscell(Qc))||(~iscell(Qz))
    error('Qc and Qz must be cell arrays.');
end
if (~isnumeric(b_axis))||(~isreal(b_axis))
    error('b_axis must be a real vector.');
end         
end

% "How many divisions?"
%
% Joseph Stalin, when told that
% Vatican was a powerful entity

