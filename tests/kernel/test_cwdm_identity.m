% Compares segmented single-substance descriptors with stock WP0 records.
% Syntax:
%
%                    result=test_cwdm_identity(record_root)
%
% Parameters:
%
%     record_root - folder containing WP0 basis_*.mat records recursively
%
% Outputs:
%
%     result      - exact descriptor comparison results
%
% WP0 stores descriptors, isotopes, chemistry, and basis inputs, but not the
% interaction graph. Thus complete and IK-0 descriptors are reconstructed
% directly; graph-restricted descriptors require capture of create outputs
% before they can be independently rebuilt. Non-sphten records have no
% spherical-tensor descriptor comparison. Unsupported records are counted
% explicitly, never labelled passes. Symmetry does not change descriptor
% rows, so its settings are omitted in this descriptor-only comparison.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_identity(record_root)

% Validate the record source
grumble(record_root);

% Announce the exact comparison
fprintf('TESTING: CWDM descriptor identity (T1)\n');
result=new_test_result('kernel/cwdm_identity','CWDM descriptor identity',...
                      'The single local descriptor equals the stock row set exactly.');
files=dir(fullfile(record_root,'**','basis_*.mat'));
compared=0; reordered=0; skipped=0;

% Rebuild every independently reconstructible single-substance record
for n=1:numel(files)
    record=load(fullfile(files(n).folder,files(n).name),'payload');
    payload=record.payload;
    if ~isfield(payload,'basis')||numel(payload.chem.parts)~=1||...
       ~strcmp(payload.bas_in.formalism,'sphten-liouv')||...
       ~ismember(payload.bas_in.approximation,{'none','IK-0'})
        skipped=skipped+1; continue;
    end
    sys.magnet=0; sys.isotopes=payload.isotopes;
    inter.chem.parts=payload.chem.parts; inter.chem.concs=1;
    options=payload.bas_in;
    fields=setdiff(fieldnames(options),{'formalism','projections','longitudinal','zero_quantum'});
    for k=1:numel(fields)
        if startsWith(fields{k},'sym_')
            options=rmfield(options,fields{k});
        elseif strcmp(fields{k},'manual')
            options.manual={logical(options.manual)};
        else
            options.(fields{k})={options.(fields{k})};
        end
    end
    actual=test_spin_system(sys,inter,options);
    descriptor=actual.bas.basis{1};
    identical=isequal(descriptor,payload.basis);
    same_rows=identical||isequal(sortrows(descriptor),sortrows(payload.basis));
    compared=compared+1; reordered=reordered+(same_rows&&~identical);
    result=test_true(result,files(n).name,same_rows,...
                     'single-substance descriptor rows equal the stock descriptor exactly');
    fprintf('T1_RECORD %s exact=%d row_set=%d\n',...
            fullfile(files(n).folder,files(n).name),identical,same_rows);
end

% Require actual records rather than an empty success
result=test_true(result,'available records',compared>0,...
                 'at least one stock descriptor was independently reconstructed');
fprintf('T1_SUMMARY compared=%d reordered=%d skipped=%d available=%d\n',...
        compared,reordered,skipped,numel(files));

end

% Consistency enforcement
function grumble(record_root)
if ~ischar(record_root)||~isfolder(record_root)
    error('record_root must name a WP0 record directory.');
end
end
