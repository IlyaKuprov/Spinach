% Tests CG phases and their SLE and PCS consumers. Syntax:
%
%                      result=test_cg_consumers()
%
% Outputs:
%
%     result - regression result with explanatory messages
%
% Wigner product identities are checked against independent rotation
% matrices, with unequal left and right projections. PCS is checked
% against the point model and non-axial distributed point sources.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cg_consumers()

% Initialise the result
result=new_test_result('kernel/cg_consumers','CG consumer conventions',...
                       'Preserve CG phases, Wigner projections, and multipole PCS.');

% Check antisymmetric coupling before two coefficients can cancel its phase
result=test_close(result,'singlet reversed projections',...
                  clebsch_gordan(0,0,0.5,-0.5,0.5,0.5),-1/sqrt(2),1e-14,0,...
                  'reversing the singlet projections reverses its sign');
for mag=-1:1
    direct=clebsch_gordan(2,mag,2,mag,1,0);
    swapped=clebsch_gordan(2,mag,1,0,2,mag);
    result=test_close(result,'unequal-spin exchange phase',direct,-swapped,1e-14,0,...
                      'odd L1+L2-L changes the sign when the coupled pairs are exchanged');
end

% Check maximal coupling against its factorial-free binomial expression
result=test_close(result,'maximal coupling binomial',clebsch_gordan(8,0,4,0,4,0),...
                  sqrt(nchoosek(8,4)^2/nchoosek(16,8)),1e-14,0,...
                  'maximal total spin agrees with the binomial coefficient formula');

% Build fresh SLE matrices without accepting an existing cached answer
cache_root=fileparts(which('sle_operators'));
cache_file=fullfile(cache_root,'sle_operators_rank_4_int_3_1_2.mat');
cache_exists=exist(cache_file,'file')==2;
if cache_exists
    backup_file=tempname;
    movefile(cache_file,backup_file);
    cleanup=onCleanup(@()movefile(backup_file,cache_file,'f'));
else
    cleanup=onCleanup(@()delete(cache_file));
end
[Lx,Ly,Lz,D,space_basis]=sle_operators(4,[3 1 2]);
result=test_close(result,'rotation commutator',full(Lx*Ly-Ly*Lx),full(1i*Lz),...
                  1e-13,1e-13,'spatial rotations obey angular momentum algebra');

% Evaluate normalised basis functions at non-axial Euler angles
angles=[0.31 0.73 -0.48;-0.62 1.17 0.29];
for point=1:size(angles,1)
    values=zeros(1,size(space_basis,1));
    rotations=cell(1,5);
    for rank=0:4
        rotations{rank+1}=wigner(rank,angles(point,1),angles(point,2),angles(point,3));
        rows=find(space_basis(:,1)==rank);
        for k=rows'
            values(k)=sqrt(2*rank+1)*rotations{rank+1}...
                      (rank+1-space_basis(k,2),rank+1-space_basis(k,3));
        end
    end

    % Restrict source ranks so no product terms fall outside the basis
    for rank=1:3
        sources=find(space_basis(:,1)<=4-rank);
        for m=1:2*rank+1
            for n=1:2*rank+1
                actual=values*D{rank}{m,n}(:,sources);
                expected=rotations{rank+1}(m,n)*values(sources);
                result=test_close(result,'Wigner product identity',actual,expected,1e-13,1e-13,...
                                  'left and right projections reproduce the full Wigner product');
            end
        end
    end
end

% Check the point-density limit with off-diagonal susceptibility
nxyz=[4 5 7;-6 4 8;9 -7 5;3 -8 -6]; mxyz=[0.1 -0.2 0.3];
chi=[0.7 0.2 -0.3;0.2 -0.1 0.4;-0.3 0.4 -0.6];
actual=lpcs(nxyz,mxyz,0,{1/sqrt(4*pi)},chi);
expected=ppcs(nxyz,mxyz,chi);
result=test_close(result,'PCS point limit',actual,expected,1e-12,1e-13,...
                  'rank-zero density reproduces point pseudocontact shifts');

% Compare a converged non-axial multipole expansion with discrete point PCS
xyz=mxyz+[0.21 -0.13 0.17;-0.16 0.23 0.11;0.09 0.07 -0.19];
weights=[0.2;0.3;0.5]; ranks=0:14;
moments=points2mult(xyz,mxyz,weights,ranks,'points');
actual=lpcs(nxyz,mxyz,ranks,moments,chi);
expected=zeros(size(actual));
for n=1:numel(weights)
    expected=expected+weights(n)*ppcs(nxyz,xyz(n,:),chi);
end
result=test_close(result,'PCS distributed point sources',actual,expected,1e-11,1e-12,...
                  'non-axial multipoles reproduce the weighted point-source PCS');

end


