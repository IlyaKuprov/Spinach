% Tests valid and rejected permutations of the shared Diels-Alder model.
% Syntax:
%
%                       result=test_dac_symmetry()
%
% Outputs:
%
%    result - full-input symmetry guards and solvent projection checks
%
% The two proposed S2 groups fail complete interaction-tensor invariance.
% The solvent S3 projector is tested with the shipped IK-2 basis against
% unprojected solvent dynamics, including its original T1/T2 relaxation.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dac_symmetry()

% Load the unchanged molecular interactions
fprintf('TESTING: Diels-Alder permutation symmetries\n');
result=new_test_result('kernel/dac_symmetry','Diels-Alder symmetry',...
                      'Full-tensor guards and solvent S3 equivalence.');
[sys,inter,bas]=dac_reaction(); sys.magnet=14.1;
s=test_spin_system(sys,inter,bas);
result=test_true(result,'solvent-only declaration',...
                 isequal(bas.sym_group,{{},{},{},{},{'S3'}})&&...
                 isequal(bas.sym_spins,{{},{},{},{},{[1 2 3]}}),...
                 'only the solvent satisfies the complete tensor contract');

% Reject each proposed S2 independently without changing physical parameters
for n=[1 4]

    % Address symmetry spins in substance-local columns
    bad=bas; bad.sym_group{n}={'S2'};
    if n==1
        bad.sym_spins{n}={[5 6]};
    else
        bad.sym_spins{n}={[7 8]};
    end
    rejected=false;
    try
        basis(s,bad);
    catch err
        rejected=contains(err.message,'Zeeman interaction tensors are not identical');
    end
    result=test_true(result,['substance ' int2str(n) ' S2 rejected'],rejected,...
                     'equal isotropic shifts do not imply identical CSA tensors');

end

% Obtain the no-symmetry reference in exactly the same descriptor basis
plain=basis(s,rmfield(bas,{'sym_group','sym_spins'}));
result=test_true(result,'unchanged descriptors',isequal(s.bas.basis,plain.bas.basis),...
                 'symmetry stores projectors without reordering basis rows');
rows=(s.bas.offsets(5)+1):s.bas.offsets(6);
P=s.bas.sym_fact(5).irr_projectors{1};
result=test_close(result,'orthonormal solvent projector',P'*P,speye(size(P,2)),...
                  1e-12,0,'SALC columns are orthonormal');
result=test_true(result,'solvent reduction',size(P,2)<numel(rows),...
                 'the fully symmetric subspace is smaller than the solvent block');

% Construct the solvent Hamiltonian from the complete production model
H=hamiltonian(assume(plain,'nmr')); H=H(rows,rows);
rho=state(plain,'L+',[28 29 30]); rho=rho(rows);
coil=coil_state(plain,'L+',[28 29 30],'exact'); coil=coil(rows);
result=test_close(result,'Hamiltonian closure',H*P,P*(P'*H*P),...
                  1e-10,1e-12,'solvent evolution preserves its symmetric subspace');
result=test_close(result,'preparation closure',P*(P'*rho),rho,...
                  1e-12,0,'collective preparation is fully symmetric');
result=test_close(result,'detection closure',P*(P'*coil),coil,...
                  1e-12,0,'collective detection is fully symmetric');

% Trace reacting spins to evaluate the original solvent relaxation cheaply
solvent=kill_spin(plain,1:27); R=relaxation(assume(solvent,'nmr'));
solvent_rows=(solvent.bas.offsets(5)+1):solvent.bas.offsets(6);
R=R(solvent_rows,solvent_rows);
result=test_close(result,'relaxation closure',R*P,P*(P'*R*P),...
                  1e-12,0,'the solvent relaxation preserves the symmetric subspace');

% Compare complex FIDs with the correct projected state and detection vector
L=H+1i*R; timing=linspace(0,1,21);
fid=zeros(size(timing)); sym_fid=zeros(size(timing));
for n=1:numel(timing)
    fid(n)=coil'*expm(full(-1i*L*timing(n)))*rho;
    sym_fid(n)=(P'*coil)'*expm(full(-1i*(P'*L*P)*timing(n)))*(P'*rho);
end
result=test_close(result,'solvent complex FID',sym_fid,fid,...
                  1e-10,1e-11,'projection preserves amplitude, phase, and decay');

end


