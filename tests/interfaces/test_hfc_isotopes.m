% Tests isotope-resolved Gaussian and ORCA hyperfine imports. Syntax:
%
%                       result=test_hfc_isotopes()
%
% Outputs:
%
%     result - regression checks for tensor scaling, provenance,
%              thresholding, purging, and unchanged NMR imports
%
% Complete tensors are compared, including anisotropic and off-diagonal
% components. The source isotopes come from the shipped electronic
% structure logs, not from assumed naturally abundant isotopes.
%
% talos@spindynamics.org

function result=test_hfc_isotopes()

% Describe the physical target
result=new_test_result('interfaces/hfc_isotopes','EPR isotope conversion',...
                       'hyperfine tensors follow the requested nuclear gyromagnetic ratio.');

% Read the Gaussian source used by the shipped ENDOR example
root_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
props=gparse(fullfile(root_dir,'examples','standard_systems','nitroxide.log'));
atom=find(strcmp(props.symbols,'N')); options.no_xyz=1;
result=test_true(result,'Gaussian source isotope',props.isotopes(atom)==14,...
                 'the printed hyperfine belongs to 14N, not the target 15N');
[sys14,inter14]=g2spinach(props,{{'E','E'},{'N','14N'}},[0 0],options);
[~,inter15]=g2spinach(props,{{'E','E'},{'N','15N'}},[0 0],options);
source_hfc=1e6*gauss2mhz(props.hfc.full.matrix{atom}/2);
result=test_close(result,'same isotope',inter14.coupling.matrix{1,2},...
                  source_hfc,0,0,'same-isotope import is unchanged');
[~,swapped]=isoswap(sys14,inter14,1,'15N');
result=test_close(result,'negative gamma ratio',inter15.coupling.matrix{1,2},...
                  swapped.coupling.matrix{1,2},1e-8,1e-14,...
                  'direct 15N import agrees with isotope replacement');
result=test_close(result,'reverse stored pair',inter15.coupling.matrix{2,1},...
                  swapped.coupling.matrix{2,1},1e-8,1e-14,...
                  'both halves of the tensor carry the isotope scaling');

% Scale Gaussian proton tensors without assuming negative gamma ratios
[~,protons]=g2spinach(props,{{'E','E'},{'H','1H'}},[0 0],options);
[~,deuterons]=g2spinach(props,{{'E','E'},{'H','2H'}},[0 0],options);
result=test_close(result,'Gaussian deuteration',cell2mat(deuterons.coupling.matrix),...
                  cell2mat(protons.coupling.matrix)*(spin('2H')/spin('1H')),1e-8,1e-14,...
                  'every proton tensor scales with the positive deuterium gamma ratio');

% Check that a larger target coupling survives the threshold and purge
options.min_hfc=norm(source_hfc,'fro')*(1+abs(spin('15N')/spin('14N')))/2;
options.purge='on';
[retained,~]=g2spinach(props,{{'E','E'},{'N','15N'}},[0 0],options);
[removed,~]=g2spinach(props,{{'E','E'},{'N','14N'}},[0 0],options);
result=test_true(result,'threshold after growth',...
                 isequal(retained.isotopes,{'15N','E'})&&isequal(removed.isotopes,{'E'}),...
                 '15N grows above the threshold while source 14N stays below it');

% Check the strict less-than threshold boundary
options.min_hfc=norm(source_hfc,'fro');
[retained,~]=g2spinach(props,{{'E','E'},{'N','14N'}},[0 0],options);
result=test_true(result,'threshold equality',numel(retained.isotopes)==2,...
                 'a tensor exactly at the existing threshold is retained');
options=rmfield(options,{'min_hfc','purge'});

% Preserve NMR conversion independently of hyperfine isotope provenance
[sys_nmr,inter_nmr]=g2spinach(props,{{'N','15N'}},0,options);
unknown=rmfield(props,'isotopes');
[sys_ref,inter_ref]=g2spinach(unknown,{{'N','15N'}},0,options);
result=test_true(result,'NMR independent of provenance',...
                 isequaln(sys_nmr,sys_ref)&&isequaln(inter_nmr,inter_ref),...
                 'NMR does not use the EPR source-isotope metadata');

% Refuse to guess the source isotope when its metadata is absent
caught=false;
try
    g2spinach(unknown,{{'E','E'},{'N','15N'}},[0 0],options);
catch exception
    caught=contains(exception.message,'not implemented')&&...
           contains(exception.message,'source isotope');
end
result=test_true(result,'missing provenance',caught,...
                 'a nonempty hyperfine tensor needs an explicit source isotope');

% Unprinted Gaussian isotope entries must not become natural-abundance guesses
unknown=props; unknown.isotopes(atom)=0; caught=false;
try
    g2spinach(unknown,{{'E','E'},{'N','15N'}},[0 0],options);
catch exception
    caught=contains(exception.message,'not implemented')&&...
           contains(exception.message,'source isotope');
end
result=test_true(result,'unprinted Gaussian isotope',caught,...
                 'zero in the parser metadata does not identify a source isotope');

% No nuclear hyperfine provenance is needed for an electron-only import
[electron,~]=g2spinach(unknown,{{'E','E'}},0,options);
result=test_true(result,'electron only',isequal(electron.isotopes,{'E'}),...
                 'the provenance requirement applies only to imported nuclear tensors');

% Read ORCA isotope strings and compare every proton tensor component
props=oparse(fullfile(root_dir,'examples','esr_liq_pulsed',...
                      'data_import','orca_methyl_radical.out'));
atoms=find(strcmp(props.symbols,'H'));
[~,protons]=g2spinach(props,{{'E','E'},{'H','1H'}},[0 0],options);
[~,deuterons]=g2spinach(props,{{'E','E'},{'H','2H'}},[0 0],options);
for n=1:numel(atoms)

    % Compare same-isotope and positive-ratio isotope conversion
    source_hfc=1e6*gauss2mhz(props.hfc.full.matrix{atoms(n)}/2);
    result=test_true(result,['ORCA isotope ' num2str(n)],...
                     strcmp(props.isotopes{atoms(n)},'1H'),...
                     'ORCA HFC provenance comes from A:ISTP, not Q:ISTP');
    result=test_close(result,['ORCA proton ' num2str(n)],...
                      protons.coupling.matrix{n,end},source_hfc,0,0,...
                      'same-isotope ORCA import is unchanged');
    result=test_close(result,['ORCA deuteron ' num2str(n)],...
                      deuterons.coupling.matrix{n,end},...
                      source_hfc*(spin('2H')/spin('1H')),1e-8,1e-14,...
                      'the full non-diagonal tensor scales by the gamma ratio');
end

% Check a shrinking target coupling before threshold and optional purge
source_norm=norm(protons.coupling.matrix{1,end},'fro');
options.min_hfc=source_norm*(1+spin('2H')/spin('1H'))/2;
[retained,pruned]=g2spinach(props,{{'E','E'},{'H','2H'}},[0 0],options);
result=test_true(result,'threshold without purge',...
                 numel(retained.isotopes)==4&&isempty(pruned.coupling.matrix{1,end}),...
                 'thresholding removes the deuterium coupling but not the nucleus');
options.purge='on';
[removed,~]=g2spinach(props,{{'E','E'},{'H','2H'}},[0 0],options);
result=test_true(result,'threshold after shrinkage',isequal(removed.isotopes,{'E'}),...
                 'all three deuterium tensors fall below the threshold');
options=rmfield(options,{'min_hfc','purge'});

% Empty ORCA isotope entries cannot supply provenance for a printed tensor
unknown=props; unknown.isotopes{atoms(1)}=[]; caught=false;
try
    g2spinach(unknown,{{'E','E'},{'H','2H'}},[0 0],options);
catch exception
    caught=contains(exception.message,'not implemented')&&...
           contains(exception.message,'source isotope');
end
result=test_true(result,'unprinted ORCA isotope',caught,...
                 'empty isotope metadata is not a licence to assume 1H');

% Partial ORCA outputs do not require metadata for unprinted tensors
props=oparse(fullfile(root_dir,'examples','visualisation','porphyrine.out'));
[partial,inter]=g2spinach(props,{{'E','E'},{'N','15N'},{'H','2H'}},[0 0 0],options);
nitrogens=find(strcmp(partial.isotopes,'15N'));
result=test_true(result,'partial ORCA output',...
                 ~isempty(nitrogens)&&all(cellfun(@isempty,inter.coupling.matrix(nitrogens,end))),...
                 'unprinted nitrogen tensors stay empty without fabricated provenance');

% Add a declared antisymmetric orbital contribution to the parsed HFC
atom=find(strcmp(props.symbols,'H'),1);
source_hfc=props.hfc.full.matrix{atom}+[0 1 -2;-1 0 3;2 -3 0]/10;
props.hfc.full.matrix{atom}=source_hfc;
[~,inter]=g2spinach(props,{{'E','E'},{'H','2H'}},[0 0],options);
result=test_close(result,'nonsymmetric isotope scaling',inter.coupling.matrix{1,end},...
                  1e6*gauss2mhz(source_hfc/2)*(spin('2H')/spin('1H')),1e-8,1e-14,...
                  'a declared antisymmetric HFC part is scaled without symmetrisation');

end

