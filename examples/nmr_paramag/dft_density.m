% Simulation of the pseudocontact shift field of the Europium(III) complex of
% 1,4,7,10-tetrakis(2-pyridylmethyl)-1,4,7,10-tetraazacyclododecane. The spin
% density, the hyperfine couplings and the susceptibility tensor are imported
% from a DFT calculation. The partial differential equation used for the delo-
% calised model solution is described in:
%
%                    http://dx.doi.org/10.1039/C4CP03106G
%
% One outlier point is due to the presence of contact shift due to the isotro-
% pic hyperfine coupling being non-sero for that particular nucleus. Point mo-
% del and Kuprov equation do not include contact shifts.
%
% g.t.p.charnock@oerc.ox.ac.uk
% ilya.kuprov@weizmann.ac.il

function dft_density()

% Load unpaired electron probability density
load('tetra_py_probden.mat','ext','xyz','probden','dx');

% Load DFT data (HFCs are read in Gauss)
props=gparse('tetra_py_dft_run.log');

% Normalize probability density
probden=probden/(trapz(trapz(trapz(probden)))*dx^3);

% Get susceptibility tensor
[~,~,rank2]=mat2sphten(props.chi);
chi=sphten2mat(0,[0 0 0],rank2);
              
% Get point pseudocontact shifts
pcs_point=ppcs(xyz,[0 0 0],chi);

% Solve Kuprov equation
[pcs_distr,pcs_3d]=kpcs(probden,chi,ext,xyz,'fft');

% Get DFT pseudocontact shifts
for n=1:(props.natoms-1) %#ok<*AGROW>
    
    % Assign isotopes
    if strcmp(props.symbols{n},'H')
        isotope='1H';
    elseif strcmp(props.symbols{n},'C')
        isotope='13C';
    elseif strcmp(props.symbols{n},'N')
        isotope='14N';
    end
    
    % Compute DFT PCS
    pcs_hfc(n)=hfc2pcs(props.hfc.full.matrix{n},chi,isotope,6); 
    
end

% Plot DFT against distributed
figure(); plot(pcs_hfc,pcs_distr,'ro');
hold on; plot([-5 5],[-5 5],'b-'); kgrid;
kxlabel('PCS from HFC tensors, ppm'); 
kylabel('PCS from Kuprov eqn, ppm');

% Plot point against distributed
figure(); plot(pcs_point,pcs_distr,'ro');
hold on; plot([-5 5],[-5 5],'b-'); kgrid;
kxlabel('PCS from point model, ppm'); 
kylabel('PCS from Kuprov eqn, ppm');

% Plot the distributed solution
figure(); pcs_3d=sqrt(abs(pcs_3d)).*sign(pcs_3d);
volplot(pcs_3d,ext); hold on; molplot(xyz,[]); 
kgrid; ktitle('PCS field');

end

