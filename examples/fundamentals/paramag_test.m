% Paramagnetic chemical shift module tests.
%
% ilya.kuprov@weizmann.ac.il

function paramag_test()

% xyz2pms against ppcs
chi=randn(3); chi=(chi+chi')/2;
nxyz=randn(1,3); sxyz=randn(1,3);
pcs_a=ppcs(nxyz,sxyz,chi);
pcs_b=trace(xyz2pms(nxyz,sxyz,chi))/3;
if abs(pcs_a-pcs_b)<1e-6
    disp('Test 1 passed.');
else
    error('Test 1 failed.');
end

% xyz2hfc + hfc2pms against xyz2pms
chi=randn(3); chi=(chi+chi')/2;
nxyz=randn(1,3); mxyz=randn(1,3);
nel=randi(10); isotope='15N';
A=xyz2hfc(mxyz,nxyz,isotope,nel);
[~,pms_tensor_a]=hfc2pms(A,chi,isotope,nel);
pms_tensor_b=xyz2pms(nxyz,mxyz,chi);
if norm(pms_tensor_a-pms_tensor_b,'fro')<1e-6
    disp('Test 2 passed.');
else
    error('Test 2 failed.');
end

end

