% Computing PCS using different models in a basic Cu(II) porphyrin complex.
% The metal is at the origin. See the "getting started" manual at
%
% http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis
%
% The paper describing the distributed PCS model used below is available at
%
%                  http://dx.doi.org/10.1039/c6cp05437d
%
% e.suturina@soton.ac.uk

function porphyrin_example_2()

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

% Cu(II) g-tensor
g_cu=diag([2.0000 2.0000 2.2000]);

% Curie susceptibility tensor
chi_cu=g2chi(g_cu,298,1/2);

% Metal position 
mxyz=[0 0 0];

% PCS calculation using the point model
point_pcs_cu=ppcs(nxyz,mxyz,chi_cu);

% Coordinates of all atoms with significant spin population
xyz=[ 0.000000000000      0.000000000000      0.000000000000
     -1.452280156504      1.452280156496      0.000000000000
      1.452280156504     -1.452280156496      0.000000000000
     -1.452280156496     -1.452280156504      0.000000000000
      1.452280156496      1.452280156504      0.000000000000];
    
% Mulliken spin populations
rho=[0.6 0.1 0.1 0.1 0.1]';

% Multipole ranks to include
L=0:14;

% Multipole moments
Ilm=points2mult(xyz,mxyz,rho,L,'points');

% PCS calculation using the distributed model
distr_pcs_cu=lpcs(nxyz,mxyz,L,Ilm,chi_cu);

% Parse the ORCA log
props=oparse('cu_porph_hfc.out');

% Extract hyperfine tensors
hfcs=props.hfc.full.matrix(26:37);

% Compute DFT PCS
dft_pcs_cu=zeros(size(point_pcs_cu));
for n=1:numel(hfcs)
    dft_pcs_cu(n)=hfc2pcs(hfcs{n},chi_cu,'1H',1);
end

% Comparison of point with distributed with DFT
disp('PCS [point, distr, dft], ppm');
disp([point_pcs_cu distr_pcs_cu dft_pcs_cu]);

end

