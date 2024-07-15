% Projection of the zero-field triplet state with user-specified
% populations of Cartesian ZFS eigenstates onto the higher-field 
% ZFS + Zeeman eigenstates. This is commonly seen in triplet DNP
% with photo-generated two-electron triplets. Syntax:
%
%    ZFS      - 3x3 ZFS tensor (Hz) in the laboratory frame
%               of reference; use zfs2mat() to get it from
%               D, E, and molecular Euler angles
%
%    pops     - a three-element vector with populations of
%               X, Y, and Z eigenstates of the ZFS tensor
%               at zero magnetic field, order: [pX pY pZ]
%
%    Z        - 3x3 Zeeman interaction tensor (Hz/Tesla) in
%               the laboratory frame of reference; use func- 
%               tions like axrh2mat() to get it from eigen-
%               values and molecular Euler angles
%
%    B        - magnetic field directed along the Z axis of
%               the laboratory frame of reference, Tesla
%
%    idx      - index of the electron triplet (use 'E3') in
%               the sys.isotopes list
%
% Outputs:
%
%    rho      - spin density matrix (Hilbert space) or state
%               vector (Liouville space)
%
% guinevere.mathies@uni-konstanz.de
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=zftrip.m>

function rho=zftrip(spin_system,ZFS,pops,Z,B,idx)

% Check consistency
grumble(spin_system,ZFS,pops,Z,B,idx);

% Get electron spin-1 operators
A=pauli(3); Ex=A.x; Ey=A.y; Ez=A.z;

% Build ZFS Hamiltonian in the laboratory frame
H_ZFS=2*pi*(ZFS(1,1)*Ex*Ex+ZFS(1,2)*Ex*Ey+ZFS(1,3)*Ex*Ez+...
            ZFS(2,1)*Ey*Ex+ZFS(2,2)*Ey*Ey+ZFS(2,3)*Ey*Ez+...
            ZFS(3,1)*Ez*Ex+ZFS(3,2)*Ez*Ey+ZFS(3,3)*Ez*Ez);

% Diagonalise ZFS Hamiltonian
[V,D]=eig(full(H_ZFS),'vector');

% Match literature convention (|zz|>|xx|>|yy|)
[~,ord]=sort(abs(D),'ascend'); V=V(:,ord);
Vz=V(:,3); Vx=V(:,2); Vy=V(:,1);

% Build zero field density matrix from populations
DM=Vx*pops(1)*Vx'+Vy*pops(2)*Vy'+Vz*pops(3)*Vz'; 

% Build Zeeman Hamiltonian in the laboratory frame
Bx=0; By=0; Bz=B;
H_Z=2*pi*(Z(1,1)*Ex*Bx+Z(1,2)*Ex*By+Z(1,3)*Ex*Bz+...
          Z(2,1)*Ey*Bx+Z(2,2)*Ey*By+Z(2,3)*Ey*Bz+...
          Z(3,1)*Ez*Bx+Z(3,2)*Ez*By+Z(3,3)*Ez*Bz);

% Diagonalise high field Hamiltonian
[V,~]=eig(full(H_Z+H_ZFS),'vector');

% Drop high-field coherences
DM=V*diag(diag(V'*DM*V))*V';

% Construct a Spinach state
T=irr_sph_ten(3); rho=complex(0);
for n=1:numel(T), T{n}=T{n}/norm(T{n},'fro')^2; end
rho=rho+trace(T{1}'*DM)*state(spin_system,'T0,0' ,idx);
rho=rho+trace(T{2}'*DM)*state(spin_system,'T1,+1',idx);
rho=rho+trace(T{3}'*DM)*state(spin_system,'T1,0' ,idx);
rho=rho+trace(T{4}'*DM)*state(spin_system,'T1,-1',idx);
rho=rho+trace(T{5}'*DM)*state(spin_system,'T2,+2',idx);
rho=rho+trace(T{6}'*DM)*state(spin_system,'T2,+1',idx);
rho=rho+trace(T{7}'*DM)*state(spin_system,'T2,0' ,idx);
rho=rho+trace(T{8}'*DM)*state(spin_system,'T2,-1',idx);
rho=rho+trace(T{9}'*DM)*state(spin_system,'T2,-2',idx);

end

% Consistency enforcement
function grumble(spin_system,ZFS,pops,Z,B,idx)
if (~isnumeric(ZFS))||(~isreal(ZFS))||(size(ZFS,1)~=3)||...
   (size(ZFS,2)~=3)||(norm(ZFS-ZFS','fro')>1e-6*norm(ZFS,'fro'))
    error('ZFS must be a real symmetric 3x3 matrx.');
end
if (~isnumeric(Z))||(~isreal(Z))||(size(Z,1)~=3)||...
   (size(Z,2)~=3)||(norm(Z-Z','fro')>1e-6*norm(Z,'fro'))
    error('Z must be a real symmetric 3x3 matrx.');
end
if (~isnumeric(pops))||(~isreal(pops))||(numel(pops)~=3)||(sum(pops)~=1)
    error('pops must be a real three-element vector with a unit sum.');
end
if (~isnumeric(B))||(~isreal(B))||(~isscalar(B))
    error('B must be a real scalar.');
end
if (~isnumeric(idx))||(~isreal(idx))||(~isscalar(idx))||...
   (idx<1)||(mod(idx,1)~=0)
    error('idx must be a real positive integer.');
end
if idx>numel(spin_system.comp.isotopes)
    error('idx exceeds the number of particles in the system.');
end
if ~strcmp(spin_system.comp.isotopes{idx},'E3')
    error(['spin ' int2str(idx) ' is not a triplet electron.']);
end
end

% Lupus dentis, taurus cornis.
%
% A Latin proverb

