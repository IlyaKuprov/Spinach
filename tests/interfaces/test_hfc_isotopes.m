% Test: interfaces/hfc_isotopes | EPR isotope conversion
% Tests Gaussian and ORCA hyperfine import and isotope substitution. Syntax:
%
%                       result=test_hfc_isotopes()
%
% Outputs:
%
%     result - regression checks for source isotope labels, tensor values,
%              thresholding, purging, and subsequent isotope substitutions
%
% talos@spindynamics.org

function result=test_hfc_isotopes()

% Describe the physical target
result=new_test_result('interfaces/hfc_isotopes','EPR isotope conversion',...
                       'hyperfine tensors belong to the isotope in the log until explicitly rescaled');

% Compare the nitrogen coupling in the shipped ENDOR example with the Gaussian tensor
root_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
props=gparse(fullfile(root_dir,'examples','standard_systems','nitroxide.log'));
atom=find(strcmp(props.symbols,'N'));
source_hfc=1e6*gauss2mhz(props.hfc.full.matrix{atom}/2);
options.no_xyz=1;
[sys,inter]=g2spinach(props,{'E','N'},[0 0],options);
result=test_true(result,'Gaussian isotope',isequal(sys.isotopes,{'14N','E'}),...
                 'Gaussian labels the nitrogen hyperfine tensor N(14)');
result=test_close(result,'Gaussian tensor',inter.coupling.matrix{1,2},source_hfc,0,0,...
                  'the parsed tensor is converted from Gauss to Hz and stored with half weight');

% A negative magnetogyric ratio reverses both stored pair tensors
[sys,inter]=isoswap(sys,inter,1,'15N');
scaled_hfc=source_hfc*(spin('15N')/spin('14N'));
result=test_true(result,'nitrogen substitution',isequal(sys.isotopes,{'15N','E'}),...
                 'only the nitrogen isotope is replaced');
result=test_close(result,'negative gamma ratio',[inter.coupling.matrix{1,2} inter.coupling.matrix{2,1}],...
                  [scaled_hfc scaled_hfc],0,1e-14,'both directions acquire the signed gamma ratio');

% A tensor on the threshold survives; raising the threshold removes only the coupling
options.min_hfc=norm(source_hfc,'fro');
[~,inter]=g2spinach(props,{'E','N'},[0 0],options);
result=test_close(result,'threshold equality',inter.coupling.matrix{1,2},source_hfc,0,0,...
                  'a tensor at the threshold is retained unchanged');
options.min_hfc=options.min_hfc+eps(options.min_hfc);
[sys,inter]=g2spinach(props,{'E','N'},[0 0],options);
result=test_true(result,'threshold without purge',isequal(sys.isotopes,{'14N','E'})&&...
                 isempty(inter.coupling.matrix{1,2})&&isempty(inter.coupling.matrix{2,1}),...
                 'both directions are removed but the nitrogen remains');

% Purging also removes the uncoupled nitrogen
options.purge='on';
[sys,~]=g2spinach(props,{'E','N'},[0 0],options);
result=test_true(result,'purge',isequal(sys.isotopes,{'E'}),'only the electron remains');

% A missing nuclear label matters only if that element is requested
props.isotopes{atom}='';
caught=false;
try
    g2spinach(props,{'E','N'},[0 0]);
catch exception
    caught=contains(exception.message,'props.isotopes')&&strcmp(exception.stack(1).name,'grumble');
end
result=test_true(result,'missing isotope',caught,'the grumbler rejects an unlabelled requested atom');
[sys,inter]=g2spinach(props,{'E'},0);
result=test_true(result,'electron only',isequal(sys.isotopes,{'E'})&&...
                 isequal(inter.zeeman.matrix,{props.g_tensor.matrix}),...
                 'the electron retains its g-tensor without a nuclear isotope');

% ORCA labels each atom, but prints hyperfine blocks only for the protons in this log
props=oparse(fullfile(root_dir,'examples','visualisation','porphyrine.out'));
[sys,inter]=g2spinach(props,{'E','Cu','N','C','H'},zeros(1,5));
protons=find(strcmp(props.symbols,'H'));
result=test_true(result,'ORCA isotopes',isequal(sys.isotopes,[props.isotopes(:)' {'E'}])&&...
                 all(cellfun(@ischar,sys.isotopes))&&all(~cellfun(@isempty,sys.isotopes))&&...
                 all(strcmp(sys.isotopes(protons),'1H'))&&all(~cellfun(@isempty,props.hfc.full.matrix(protons))),...
                 'printed A:ISTP=1 identifies the protons; all other centres have isotope labels');
source_hfc=1e6*gauss2mhz(cell2mat(props.hfc.full.matrix(protons))/2);
result=test_close(result,'ORCA proton tensors',cell2mat(inter.coupling.matrix(protons,end)'),...
                  source_hfc,0,0,'the full proton hyperfine tensors are imported as printed');
result=test_true(result,'ORCA unprinted tensors',all(cellfun(@isempty,inter.coupling.matrix(~strcmp(props.symbols,'H'),end))),...
                 'assigning an isotope does not invent a hyperfine tensor');

% A positive magnetogyric ratio reduces the proton tensors on deuteration
[sys,inter]=isoswap(sys,inter,protons,'2H');
result=test_true(result,'deuteration',all(strcmp(sys.isotopes(protons),'2H')),...
                 'the selected protons are replaced by deuterons');
result=test_close(result,'positive gamma ratio',cell2mat(inter.coupling.matrix(protons,end)'),...
                  source_hfc*(spin('2H')/spin('1H')),0,1e-14,'deuteration scales every tensor component');

end


