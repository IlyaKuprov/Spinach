% Adaptively recursed Voitlander integrator. Computes an approxi-
% mation of an integral of field-swept EPR transition over a sph-
% erical triangle. Syntax:
%
% spec=voitlander(spin_system,parameters,r1,r2,r3,tf1,tf2,tf3,...
%                 tm1,tm2,tm3,tw1,tw2,tw3,Ic,Iz,Qc,Qz,Hmw,b_axis)
%
% Parameters:
%
%     r1,r2,r3    - Cartesian coordinates of the corners of the
%                   spherical triangle, unit column vectors
%
%     tf1,tf2,tf3 - transition frequencies at the corners of the
%                   spherical triangle, real column vectors, one
%                   element per transition
%
%     tm1,tm2,tm3 - transition moments at the corners of the sphe-
%                   rical triangle, positive column vectors, one
%                   element per transition
%
%     tw1,tw2,tw3 - transition widths at the corners of the sphe-
%                   rical triangle, positive column vectors, one
%                   element per transition
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
% Outputs:
%
%     spec        - ESR spectrum integral over the triangle, array
%                   of the same dimension as b_axis
%               
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=voitlander.m>

function spec=voitlander(spin_system,parameters,r1,r2,r3,tf1,tf2,tf3,...
                         tm1,tm2,tm3,tw1,tw2,tw3,Ic,Iz,Qc,Qz,Hmw,b_axis)
                     
% Check consistency
grumble(r1,r2,r3,tf1,tf2,tf3,tm1,tm2,tm3,...
        tw1,tw2,tw3,Ic,Iz,Qc,Qz,b_axis);

% Subdivide the triangle
[r12,r23,r31]=sphtrsubd(r1,r2,r3);

% Characterise new vertices
[phi,elev]=cart2sph(r12(1),r12(2),r12(3));
theta=pi/2-elev; parameters.orientation=[0 theta phi];
[tf12,tm12,tw12]=eigenfields(spin_system,parameters,Iz,Qz,Ic,Qc,Hmw);

[phi,elev]=cart2sph(r23(1),r23(2),r23(3)); 
theta=pi/2-elev; parameters.orientation=[0 theta phi];
[tf23,tm23,tw23]=eigenfields(spin_system,parameters,Iz,Qz,Ic,Qc,Hmw);

[phi,elev]=cart2sph(r31(1),r31(2),r31(3));
theta=pi/2-elev; parameters.orientation=[0 theta phi];
[tf31,tm31,tw31]=eigenfields(spin_system,parameters,Iz,Qz,Ic,Qc,Hmw);

% Compute the subdivided integral
spec_sub=trint(tf1,tf12,tf31,tm1,tm12,tm31,tw1,tw12,tw31,r1,r12,r31,b_axis)+...
         trint(tf12,tf2,tf23,tm12,tm2,tm23,tw12,tw2,tw23,r12,r2,r23,b_axis)+...
         trint(tf31,tf23,tf3,tm31,tm23,tm3,tw31,tw23,tw3,r31,r23,r3,b_axis)+...
         trint(tf12,tf23,tf31,tm12,tm23,tm31,tw12,tw23,tw31,r12,r23,r31,b_axis);
     
% Compute the direct integral
spec_dir=trint(tf1,tf2,tf3,tm1,tm2,tm3,tw1,tw2,tw3,r1,r2,r3,b_axis);

% If the accuracy is insufficient, make a recursive call
if norm(spec_dir-spec_sub,inf)>parameters.int_tol
    spec=voitlander(spin_system,parameters,r1,r12,r31,tf1,tf12,tf31,tm1,tm12,tm31,tw1,tw12,tw31,Ic,Iz,Qc,Qz,Hmw,b_axis)+...
         voitlander(spin_system,parameters,r12,r2,r23,tf12,tf2,tf23,tm12,tm2,tm23,tw12,tw2,tw23,Ic,Iz,Qc,Qz,Hmw,b_axis)+...
         voitlander(spin_system,parameters,r31,r23,r3,tf31,tf23,tf3,tm31,tm23,tm3,tw31,tw23,tw3,Ic,Iz,Qc,Qz,Hmw,b_axis)+...
         voitlander(spin_system,parameters,r12,r23,r31,tf12,tf23,tf31,tm12,tm23,tm31,tw12,tw23,tw31,Ic,Iz,Qc,Qz,Hmw,b_axis);
else
    spec=spec_sub;
end

end

% Single triangle integrator
function spec=trint(tf1,tf2,tf3,tm1,tm2,tm3,...
                    tw1,tw2,tw3, r1, r2, r3,b_axis)
                 
% Preallocate the spectrum
spec=zeros(size(b_axis),'like',1i);

% Area of spherical triangle
S=sphtarea(r1,r2,r3);

% Ignore hanging transitions
ntrans=min([numel(tf1) numel(tf2) numel(tf3)]);
if numel(tf1)>ntrans
    [~,kill_mask]=mink(tm1,numel(tm1)-ntrans); tf1(kill_mask)=[]; tm1(kill_mask)=[];
end
if numel(tf2)>ntrans
    [~,kill_mask]=mink(tm2,numel(tm2)-ntrans); tf2(kill_mask)=[]; tm2(kill_mask)=[];
end
if numel(tf3)>ntrans
    [~,kill_mask]=mink(tm3,numel(tm3)-ntrans); tf3(kill_mask)=[]; tm3(kill_mask)=[];
end

% Build spectrum
for k=1:ntrans

    % Get signal information
    min_freq=min([tf1(k) tf2(k) tf3(k)]);
    max_freq=max([tf1(k) tf2(k) tf3(k)]);
    tran_mom=mean([tm1(k) tm2(k) tm3(k)]);
    line_width=mean([tw1(k) tw2(k) tw3(k)]);

    % Find the relevant part of the axis
    b_mask=(b_axis>min_freq-3*line_width)&(b_axis<max_freq+3*line_width);

    % Convolutions of Lorentzians with triangles
    spec(b_mask)=spec(b_mask)+S*lorentzcon([tf1(k) tf2(k) tf3(k)],...
                                tran_mom,line_width,b_axis(b_mask));
                  
end

end

% Consistency enforcement
function grumble(r1,r2,r3,tf1,tf2,tf3,tm1,tm2,tm3,...
                 tw1,tw2,tw3,Ic,Iz,Qc,Qz,b_axis)
if (~isnumeric(r1))||(~isreal(r1))||(~iscolumn(r1))||...
   (numel(r1)~=3)||(abs(norm(r1,2)-1)>1e-6)||...
   (~isnumeric(r2))||(~isreal(r2))||(~iscolumn(r2))||...
   (numel(r2)~=3)||(abs(norm(r2,2)-1)>1e-6)||...
   (~isnumeric(r3))||(~isreal(r3))||(~iscolumn(r3))||...
   (numel(r3)~=3)||(abs(norm(r3,2)-1)>1e-6)
    error('r1,r2,r3 must be real unit column 3-vectors.');
end
if (~isnumeric(tf1))||(~isreal(tf1))||(~iscolumn(tf1))||...
   (~isnumeric(tf2))||(~isreal(tf2))||(~iscolumn(tf2))||...
   (~isnumeric(tf3))||(~isreal(tf3))||(~iscolumn(tf3))
    error('tf1,tf2,tf3 must be real column vectors.');
end
if (~isnumeric(tm1))||(~isreal(tm1))||(~iscolumn(tm1))||...
   (~isnumeric(tm2))||(~isreal(tm2))||(~iscolumn(tm2))||...
   (~isnumeric(tm3))||(~isreal(tm3))||(~iscolumn(tm3))
    error('tm1,tm2,tm3 must be real column vectors.');
end
if (~isnumeric(tw1))||(~isreal(tw1))||(~iscolumn(tw1))||...
   (~isnumeric(tw2))||(~isreal(tw2))||(~iscolumn(tw2))||...
   (~isnumeric(tw3))||(~isreal(tw3))||(~iscolumn(tw3))
    error('tw1,tw2,tw3 must be real column vectors.');
end
if (~isnumeric(Ic))||(size(Ic,1)~=size(Ic,2))||...
   (~isnumeric(Iz))||(size(Iz,1)~=size(Iz,2))
    error('Ic and Iz must be square Hermitian matrices.');
end
if (~iscell(Qc))||(~iscell(Qz))
    error('Qc and Qz must be cell arrays.');
end
if (~isnumeric(b_axis))||(~isreal(b_axis))
    error('b_axis must be a real vactor.');
end         
end

% "How many divisions?"
%
% Joseph Stalin, when told that
% Vatican was a powerful entity

