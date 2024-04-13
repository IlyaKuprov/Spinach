% Computes PCS using different models in basic Cu(II) porphyrin complex. See 
% the "getting started" manual at
%
% http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis
%
% The paper describing the distributed PCS model used below is available at
%
%                   http://dx.doi.org/10.1039/c6cp05437d
%
% Calculation time: minutes, 64GB of RAM required
%
% e.suturina@soton.ac.uk

function porphyrin_example_3()

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

% Cu(II) g-tensor eigenvalues
g_cu=diag([2.0000 2.0000 2.2000]);

% Curie susceptibility tensor
chi_cu=g2chi(g_cu,298,1/2);

% Parse ORCA log
props=oparse('cu_porph_hfc.out');

% Extract hyperfine tensors
hfcs=props.hfc.full.matrix(26:37);

% Compute HFC PCS
hfc_pcs_cu=zeros(numel(hfcs),1);
for n=1:numel(hfcs)
    hfc_pcs_cu(n)=hfc2pcs(hfcs{n},chi_cu,'1H',1);
end

% Parse ORCA cube and pad the density with zeros to avoid PBC effects 
pad_size=2; [density,ext]=ocparse('cu_porph_sd120.spindens.3d',pad_size);

% Compute extentens for original density without padding
zoom_vol=[pad_size/(2*pad_size+1) (pad_size+1)/(2*pad_size+1)...
          pad_size/(2*pad_size+1) (pad_size+1)/(2*pad_size+1)...
          pad_size/(2*pad_size+1) (pad_size+1)/(2*pad_size+1)];

% Draw the probability density schematic
[zoom_den,zoom_ext]=zoom_3d(density,ext,zoom_vol);
figure(); volplot(zoom_den,zoom_ext,[0.05 0.05]);
hold on; molplot(props.std_geom,[]); drawnow();

% Solve Kuprov equation
[pde_pcs_cu,pcs_cube]=kpcs(density,chi_cu,ext,nxyz,'fft');

% Compare HFC PCS with PDE PCS
disp('Pseudocontact shifts [hfc, pde], ppm');
disp([hfc_pcs_cu pde_pcs_cu]);

% Plot the PCS field schematic
[zoom_pcs,zoom_ext]=zoom_3d(pcs_cube,ext,zoom_vol);
figure(); volplot(zoom_pcs,zoom_ext,[0.02 0.02]); 
hold on; molplot(props.std_geom,[]);

end

