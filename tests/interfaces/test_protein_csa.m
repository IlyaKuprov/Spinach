% Tests canonical protein CSA import and peptide N/H overrides. Syntax:
%
%                       result=test_protein_csa()
%
% Outputs:
%
%     result - regression checks for full non-symmetric tensors,
%              PDB serial matching, BMRB shifts, selection, and validation
%
% The two-residue fixture has non-consecutive serials, a removed oxygen,
% repeated CA labels, and one unassigned proton. Each tensor is traceless
% and non-symmetric, so transposition and symmetrisation cannot pass.
%
% talos@spindynamics.org

function result=test_protein_csa()

% Describe the import contract
result=new_test_result('interfaces/protein_csa','Protein CSA file import',...
                       'imported anisotropy preserves Gaussian rows and BMRB isotropic shifts.');

% Create disposable input files
work_dir=tempname; mkdir(work_dir);
cleanup=onCleanup(@()rmdir(work_dir,'s'));
pdb_file=fullfile(work_dir,'fixture.pdb');
bmrb_file=fullfile(work_dir,'fixture.bmrb');
csa_file=fullfile(work_dir,'fixture.txt');

% Specify a non-collinear peptide bond with an unassigned proton
serials=[10 20 25 30 40 50 60];
labels={'CA','C','O','N','H','CA','HA2'};
resnums=[1 1 1 2 2 2 2];
coords=[-1 0 0;0 0 0;0 -1 0;0 1.3 0;0 1.3 1;1.4 1.3 0;1.4 2.3 0];
pdb_lines=strings(7,1); bmrb_lines=strings(6,1);
for n=1:7

    % Write records in the formats accepted by the production readers
    pdb_lines(n)=sprintf('ATOM %d %s GLY A %d %.3f %.3f %.3f 1.00 0.00',...
                         serials(n),labels{n},resnums(n),coords(n,:));
    if n<=6
        bmrb_lines(n)=sprintf('%d %d GLY %s X %.4f 0 1',...
                              n,resnums(n),labels{n},n+0.125);
    end
end
writelines(pdb_lines,pdb_file); writelines(bmrb_lines,bmrb_file);

% Supply distinct tensors in reverse order to rule out positional mapping
csa_lines=repmat("# canonical CSA fixture",33,1);
tensor=[1 2 3;4 -3 5;6 7 2];
for n=1:7

    % Keep each full Gaussian row exactly as printed
    atom=8-n;
    csa_lines(4*n+2)=sprintf('%d %s GLY %d',serials(atom),labels{atom},resnums(atom));
    for k=1:3
        csa_lines(4*n+2+k)=sprintf('%12.4f %12.4f %12.4f',atom*tensor(k,:));
    end
end
writelines(csa_lines,csa_file);

% Compare the legacy path directly with the existing CSA estimator
options.select='all'; options.pdb_mol=1; options.noshift='keep';
options.deuterate={};
output=evalc('[sys,legacy]=protein(pdb_file,bmrb_file,options);'); %#ok<NASGU>
retained=[1 2 4 5 6 7];
[~,guessed]=evalc('guess_csa_pro(resnums(retained)'',labels(retained)'',num2cell(coords(retained,:),2),options)');
result=test_true(result,'legacy guessing',isequaln(legacy.zeeman.matrix,guessed'),...
                 'without csa_file, the existing estimator still supplies the tensors');

% Import every retained carbon, nitrogen, and proton without guessing
options.csa_file=csa_file;
output=evalc('[imported,inter]=protein(pdb_file,bmrb_file,options);'); %#ok<NASGU>
result=test_true(result,'unchanged system',isequal(sys,imported),...
                 'CSA import does not change spin labels or isotopes');
result=test_close(result,'BMRB and missing shifts',cell2mat(inter.zeeman.scalar),...
                  [1.125 2.125 4.125 5.125 6.125 0],0,0,...
                  'assigned shifts come from BMRB and the unassigned spin keeps its legacy shift');
for n=1:numel(retained)

    % Compare all nine components including the antisymmetric part
    result=test_close(result,['Gaussian rows ' num2str(n)],inter.zeeman.matrix{n},...
                      retained(n)*tensor,0,0,'serial and all labels identify the tensor, not table order');
end

% Retain BMRB isotropy and imported carbons under each explicit N/H override
for mode={'tcb','bax','pol'}

    % Capture the warning without suppressing or leaking it into the runner log
    options.nh_csa=mode{1}; lastwarn('');
    output=evalc('[~,overridden]=protein(pdb_file,bmrb_file,options);'); %#ok<NASGU>
    [~,warn_id]=lastwarn;
    [~,guessed]=evalc('guess_csa_pro(resnums(retained)'',labels(retained)'',num2cell(coords(retained,:),2),options)');
    result=test_true(result,['override warning ' mode{1}],strcmp(warn_id,'protein:nh_csa_override'),...
                     'explicit nh_csa warns that imported amide N/H anisotropy is overridden');
    result=test_close(result,['override N ' mode{1}],overridden.zeeman.matrix{3},guessed{3},0,0,...
                      'amide nitrogen uses the existing estimator');
    result=test_close(result,['override H ' mode{1}],overridden.zeeman.matrix{4},guessed{4},0,0,...
                      'amide proton uses the existing estimator');
    result=test_true(result,['override remainder ' mode{1}],...
                     isequal(overridden.zeeman.matrix([1 2 5 6]),inter.zeeman.matrix([1 2 5 6]))&&...
                     isequal(overridden.zeeman.scalar,inter.zeeman.scalar),...
                     'carbon and non-amide proton tensors and all isotropic shifts are unchanged');
end
options=rmfield(options,'nh_csa');

% Select PDB serials rather than post-filter indices and deuterate after matching
options.select=[40 50]; options.deuterate={'H'};
output=evalc('[selected,selected_inter]=protein(pdb_file,bmrb_file,options);'); %#ok<NASGU>
result=test_true(result,'selection and deuteration',...
                 isequal(selected.labels,{'D','CA'})&&isequal(selected.isotopes,{'2H','13C'})&&...
                 isequal(selected_inter.zeeman.matrix,inter.zeeman.matrix([4 5])),...
                 'matching precedes isotope replacement and selection uses the original PDB serials');

% Missing entries for deleted or unselected atoms do not block the import
options.select='all'; options.deuterate={}; options.noshift='delete';
writelines(csa_lines([1:5 10:end]),csa_file);
output=evalc('[selected,selected_inter]=protein(pdb_file,bmrb_file,options);'); %#ok<NASGU>
result=test_true(result,'unassigned deletion',numel(selected.labels)==5&&...
                 isequal(selected_inter.zeeman.matrix,inter.zeeman.matrix(1:5)),...
                 'an unassigned deleted atom does not require a CSA record');
options.select=[10 20]; options.noshift='keep';
writelines(csa_lines([1:5 26:end]),csa_file);
output=evalc('[~,selected_inter]=protein(pdb_file,bmrb_file,options);'); %#ok<NASGU>
result=test_true(result,'selected-only table',...
                 isequal(selected_inter.zeeman.matrix,inter.zeeman.matrix(1:2)),...
                 'only retained PDB atoms require matching records');

% Reject malformed blocks, ambiguous identifiers, and incomplete retained data
options.select='all';
invalids=repmat({csa_lines},1,13);
invalids{1}(1)="not a comment";
invalids{2}(end)=[];
invalids{3}(7)="1 2";
invalids{4}(7)="1 2 3 4";
invalids{5}(7)="NaN 2 3";
invalids{6}(7)="Inf 2 3";
invalids{7}(6)="50 CA GLY 2";
invalids{8}=csa_lines([1:5 10:end]);
invalids{9}(6)="60 HA3 GLY 2";
invalids{10}(7)="8 14 21";
invalids{11}(6)="60 HA2 GLY";
invalids{12}(6)="60.5 HA2 GLY 2";
invalids{13}(7)="7i 14 21";
messages={'five #','four-line','three finite','three finite','three finite',...
          'three finite','duplicate','missing','labels do not match','traceless',...
          'atom header','finite integers','three finite'};
for n=1:numel(invalids)

    % Require an informative production error for each malformed input
    writelines(invalids{n},csa_file); caught=false;
    try
        output=evalc('protein(pdb_file,bmrb_file,options);'); %#ok<NASGU>
    catch exception
        caught=contains(exception.message,messages{n});
    end
    result=test_true(result,['invalid table ' num2str(n)],caught,...
                     'invalid or missing CSA data cannot silently become guessed tensors');
end

% Reject invalid file options at the input boundary
for invalid={17,"filename",'',fullfile(work_dir,'absent.txt')}

    % Exercise the character row-vector and file existence contract
    options.csa_file=invalid{1}; caught=false;
    output=evalc(['try; protein(pdb_file,bmrb_file,options); catch exception; '...
                  'caught=contains(exception.message,''options.csa_file''); end']); %#ok<NASGU>
    result=test_true(result,'invalid file option',caught,...
                     'csa_file must be a non-empty character row vector naming an existing file');
end

end


