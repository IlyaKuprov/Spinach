% Computing PCS using different models in basic Cu(II) and Co(II) porphyrin 
% complexes. The metal is at the origin. See the "getting started" manual at
%
% http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis
%
% e.suturina@soton.ac.uk

function porphyrin_example_1()

% Porphyrin ring proton coordinates
nxyz=[ 4.551635888      2.658552774      0.000000000
       2.658552774      4.551635889      0.000000000
      -2.658552774      4.551635888      0.000000000
      -4.551635889      2.658552774      0.000000000
      -4.551635888     -2.658552774      0.000000000
      -2.658552774     -4.551635889      0.000000000
       2.658552774     -4.551635888      0.000000000
       4.551635889     -2.658552774      0.000000000
       4.533874147      0.000000000      0.000000000
       0.000000000     -4.533874147      0.000000000
      -4.533874147      0.000000000      0.000000000
       0.000000000      4.533874147      0.000000000];

% Co(II) g-tensor
g_co=diag([3.0 3.0 2.0]);

% Cu(II) g-tensor
g_cu=diag([2.0 2.0 2.2]);

% Curie susceptibility tensors
chi_co=g2chi(g_co,298,1/2);
chi_cu=g2chi(g_cu,298,1/2);

% Metal position 
mxyz=[0 0 0];

% PCS calculation
point_pcs_co=ppcs(nxyz,mxyz,chi_co);
point_pcs_cu=ppcs(nxyz,mxyz,chi_cu);

% Output
disp('PCS [Co, Cu], ppm');
disp([point_pcs_co point_pcs_cu]);

end

