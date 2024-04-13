% Test of matrix exponential differentiation of second order
% Magnus product quadrature (trapdiff.m) with the result com-
% pared to the central finite difference derivative. General
% coherent + non-symmetric dissipative case is tested.
%
% u.rasulov@soton.ac.uk
% i.kuprov@soton.ac.uk

function dirdiff_2()

% Bootstrap Spinach
spin_system=bootstrap();

% Left and right drift generators, dissipative
Hd={randn(50)+1i*randn(50), randn(50)+1i*randn(50)};

% Control operator
Hc=randn(50)+1i*randn(50);

% A reasonable time step estimate
dt=mean([1/norm(Hd{1},2), 1/norm(Hd{2},2)]);

% Reasonable controls
cL=0.1*randn()*norm(Hd{1},2);
cR=0.1*randn()*norm(Hd{2},2);

% Get analytical derivatives
[DL_anl,DR_anl]=trapdiff(spin_system,Hd,Hc,dt,cL,cR);

% Get numerical derivatives
dc=sqrt(eps('double'));
H_dir_L=(1/2)*Hc+1i*dt*(sqrt(3)/12)*(Hc*Hd{2}-Hd{2}*Hc);
H_dir_R=(1/2)*Hc+1i*dt*(sqrt(3)/12)*(Hd{1}*Hc-Hc*Hd{1});
H=(Hd{1}+Hd{2})/2+cL*H_dir_L+cR*H_dir_R;
DL_num=(expm(-1i*dt*(H+dc*H_dir_L))-...
        expm(-1i*dt*(H-dc*H_dir_L)))/(2*dc);
DR_num=(expm(-1i*dt*(H+dc*H_dir_R))-...
        expm(-1i*dt*(H-dc*H_dir_R)))/(2*dc);

% Test analytical against finite difference
if (norm(DL_anl-DL_num,2)<10*sqrt(eps('double')))&&...
   (norm(DR_anl-DR_num,2)<10*sqrt(eps('double')))
   disp('test passed');
else
   error('test failed');
end

end

