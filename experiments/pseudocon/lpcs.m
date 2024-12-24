% Computes PCS from the multipole moments of the paramagnetic 
% centre probability density as described in 
%
%             http://dx.doi.org/10.1039/c6cp05437d
%
% Equation 33 is used under the assumption that all nuclei are
% outside the bounding sphere shown in Figure 1 of the paper.
%
% Syntax:
%
%            theo_pcs=lpcs(nxyz,mxyz,ranks,Ilm,chi)
% 
% Parameters: 
%
%     nxyz  - nuclear coordinates as [x y z] with multiple
%             rows, at which PCS is to be evaluated, in 
%             Angstroms.
%
%     mxyz  - paramagnetic centre coordinates as [x y z], 
%             in Angstroms.
%
%     ranks - array of multipole ranks supplied in Ilm
%
%     Ilm   - {[],[]} array of numbers corresponding to 
%             the integrals
%
%             Int[rho(r,theta,phi)*Y_lm(theta,phi)*r^l*d^3r]
%
%             for L=0 Ilm=N/2/sqrt(pi)
%             for L=1 Ilm=[imag(I11) I10 real(I11)]
%             for L=2 Ilm=[imag(I22) imag(I21) I20 ...
%                          real(I21) real(I22)]
%             
%             et cetera.
%
%     chi   - the 3x3 matrix or the five independent elements
%             of the susceptibility tensor in cubic Angstroms
%
% Output:
% 
%     theo_pcs  - predicted pseudocontact shift (in ppm) at 
%                 each of the nuclei.
%
% ilya.kuprov@weizmann.ac.uk
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=lpcs.m>

function theo_pcs=lpcs(nxyz,mxyz,ranks,Ilm,chi)

% Check consistency
grumble(nxyz,mxyz,ranks,Ilm,chi);

% Put the metal at the origin
nxyz=nxyz-kron(mxyz,ones(size(nxyz,1),1));

% Form the susceptibility tensor
if numel(chi)==5
    chi=[chi(1)     chi(2)     chi(3); 
         chi(2)     chi(4)     chi(5);
         chi(3)     chi(5)  -chi(1)-chi(4)];
end

% Compute the irreducible components
[~,~,chi2]=qform2sph(chi);

% Build look-up tables
A=zeros(max(ranks)+1,1);
for n=0:max(ranks)
    A(n+1)=cg_fast(n+2,0,n,0,2,0);
end
B=zeros(max(ranks)+1,2*(max(ranks)+1)+1,5);
for n=0:max(ranks)
    for p=-n:n
        for q=-2:2
            B(n+1,p+n+1,q+3)=cg_fast(n+2,p+q,n,p,2,q);
        end
    end
end

% Convert coordinates to spherical
[r,theta,phi]=xyz2sph(nxyz(:,1),nxyz(:,2),nxyz(:,3));

% Preallocate the answer
theo_pcs=zeros(numel(r),1);

% Use Equation 33 from http://dx.doi.org/10.1039/c6cp05437d
for i=1:numel(ranks)
    l=ranks(i);
    for m=-l:l
        for mp=-2:2
            theo_pcs=theo_pcs+sqrt(20*(2*l+1)/(pi*(2*l+5)))/12*(2*l+3)*A(l+1)*B(l+1,m+l+1,mp+3)*...
                              (sign(m)^m*Ilm{i}(abs(m)+l+1)+sqrt(-1)*(sign(m))^(m+1)*Ilm{i}(-abs(m)+l+1))*...
                              chi2(-mp+3)*r.^(-l-3).*spher_harmon(l+2,m+mp,theta,phi);
        end
    end
end

% Clean up and convert to ppm
theo_pcs=1e6*real(theo_pcs);

end

% Consistency enforcement
function grumble(nxyz,mxyz,L,Ilm,chi)
if (~isnumeric(nxyz))||(~isreal(nxyz))||(size(nxyz,2)~=3)
    error('nxyz must be an Nx3 array of atomic coordinates.');
end
if (~isnumeric(mxyz))||(size(mxyz,2)~=3)||(size(mxyz,1)~=1)||(~isreal(mxyz))
    error('mxyz parameter should be a real row vector with three elements.');
end
if (~isnumeric(L))||(~isreal(L))
    error('L must be a real vector.');
end
if (~isnumeric(chi))||(~isreal(chi))
    error('chi should be a real');
end
if (~iscell(Ilm))||(numel(Ilm)~=numel(L))
    error('Ilm should be a vell array of the same size as L')
end
for n=1:numel(L)
    l=L(n);
    if (size(Ilm{n})~=(2*l+1))
        error('Ilm{n} should have 2*L(n)+1 elements')
    end
end
end

% I am not only a rationalist, but also an empiricist. Rational thinking
% can start from false premises and be effectively superstition. But some-
% thing that was tested and found working, is working.
%
% Unknown source

