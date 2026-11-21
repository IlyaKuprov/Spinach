% Test: interfaces/s2spinach | literal spin-system reader
% Tests literal spin-system input against the bundled examples,
% independently evaluated rank-two projections, and native Hamiltonians.
% Syntax:
%
%                       result=test_s2spinach()
%
% Outputs:
%
%     result  - regression test result with explanatory messages

function result=test_s2spinach()

% State the scope of the regression
result=new_test_result('interfaces/s2spinach','literal spin-system reader',...
                       'literal records produce native inputs with matching units and axes.');

% Find the bundled decks and set up a native basis
root_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
deck_names={'ccacb.spinsys','ncac.spinsys','ncca.spinsys'};
isotope_lists={{'13C','13C','13C'},{'15N','13C','13C'},{'15N','13C','13C'}};
bas.formalism='zeeman-hilb';
bas.approximation={'none'};
for n=1:numel(deck_names)

    % Build each bundled spin system through the public interface
    deck=dir(fullfile(root_dir,'examples','standard_systems','**',deck_names{n}));
    [sys,inter]=s2spinach(fullfile(deck.folder,deck.name));
    sys.magnet=9.4;
    sys.parallel={'processes',1};
    sys.output='hush';
    spin_system=basis(create(sys,inter),bas);
    result=test_true(result,deck_names{n},...
                     isequal(spin_system.comp.isotopes,isotope_lists{n})&&...
                     all(~cellfun(@isempty,inter.zeeman.matrix))&&...
                     nnz(~cellfun(@isempty,inter.coupling.matrix))==3&&...
                     nnz(~cellfun(@isempty,inter.coupling.scalar))==2,...
                     'nuclei, shifts, and all three interacting pairs reach a native build.');

end

% Write a literal fixture with nonzero unequal angles and out-of-order records
file_name=[tempname() '.spinsys'];
file_id=fopen(file_name,'w');
remove_file=onCleanup(@()delete(file_name));
fprintf(file_id,'%s\n','spinsys {','shift 1 +8p 3.0e1p 0.3 19 37 53 # ppm',...
                       'quadrupole 1 2 1200 0.3 19 37 53','channels 2H',...
                       'nuclei 2H','}');
fclose(file_id);
[sys,inter]=s2spinach(string(file_name));

% Six independent field projections of the reference rank-two rotation
directions=[0 0 1; 1 0 0; 0 1 0; 1 1 0; 1 0 1; 0 1 1];
directions=directions./vecnorm(directions,2,2);
shift_ref=[20.4175280113463; 2.47235496505877; 1.11011702359496;...
           -3.86063743217297; -1.26546846008347; 30.4015828589815];
shift_proj=sum((directions*inter.zeeman.matrix{1}).*directions,2);
result=test_close(result,'tilted shift projections',shift_proj,shift_ref,1e-11,0,...
                  'six field directions fix the signs of all off-diagonal elements.');
quad_proj=sum((directions*inter.coupling.matrix{1}).*directions,2);
result=test_close(result,'tilted quadrupole projections',quad_proj,20*(shift_ref-8),1e-10,0,...
                  'spin one has traceless zz equal to Cq/2, with the same PAS orientation.');

% Build the spin-one coupling Hamiltonian and check its unrotated spectrum
sys.magnet=9.4;
sys.parallel={'processes',1};
sys.output='hush';
spin_system=basis(create(sys,inter),bas);
[H,Q]=hamiltonian(assume(spin_system,'labframe','couplings'));
levels=eig(full(H+orientation(Q,[0 0 0])))/(2*pi);
result=test_close(result,'quadrupole energy levels',sort(levels),[-600; 210; 390],1e-8,0,...
                  'a tilted full spin-one quadrupole preserves the principal-axis eigenvalues.');

% Use geometry to set the signed dipole in Hz without relying on the reader
[dipole,~,~,~,dip_tensor]=xyz2dd([0 0 0],[0 0 1.5],'13C','15N');
dipole_line=sprintf('dipole 2 1 %.17g 0 0 0',-dipole/(2*pi));
file_id=fopen(file_name,'w');
fprintf(file_id,'%s\n','spinsys {','channels 13C','nuclei 13C 15N',...
                       'shift 1 10p 0 0 0 0 0','shift 2 10p 0 0 0 0 0',...
                       dipole_line,'jcoupling 1 2 127 19 0 0 0 0','}');
fclose(file_id);
[sys,inter]=s2spinach(file_name);
result=test_close(result,'dipole plus anisotropic J',inter.coupling.matrix{1,2},...
                  dip_tensor/(2*pi)+diag([-19 -19 38]),1e-10,0,...
                  'the geometry fixes the signed dipole and its radians-to-cycles conversion.');

% Verify the heteronuclear conditional splitting and its Hamiltonian sign
sys.magnet=9.4;
sys.parallel={'processes',1};
sys.output='hush';
spin_system=basis(create(sys,inter),bas);
[H,Q]=hamiltonian(assume(spin_system,'nmr','couplings'));
levels=diag(full(H+orientation(Q,[0 0 0])))/(2*pi);
pair_hz=127+38+dip_tensor(3,3)/(2*pi);
result=test_close(result,'heteronuclear coupling energies',levels,...
                  pair_hz*[1;-1;-1;1]/4,1e-8,0,...
                  'one signed pair gives the expected conditional frequency splitting, not twice that value.');

% Verify the ppm sign for nuclei with opposite signed gyromagnetic ratios
H=hamiltonian(assume(spin_system,'nmr','zeeman'));
shift_hz=-sys.magnet*[spin('13C') spin('15N')]*10e-6/(2*pi);
levels=diag(full(H))/(2*pi);
result=test_close(result,'signed isotope shifts',levels,...
                  [1 1; 1 -1; -1 1; -1 -1]*shift_hz.'/2,1e-8,0,...
                  'identical ppm shifts give opposite signed frequencies.');

% Reject malformed and unsupported records at their original file line
bad_records={'dipole 1 3 10 0 0 0','s2spinach:index';...
             'shift 1 1000 5p 0 0 0 0','s2spinach:units';...
             'mixing 1 2','s2spinach:unsupported';...
             'shift 1 [expr 3] 5p 0 0 0 0','s2spinach:syntax'};
for n=1:size(bad_records,1)

    % Replace only the interaction under test
    file_id=fopen(file_name,'w');
    fprintf(file_id,'%s\n','spinsys {','nuclei 13C 15N',bad_records{n,1},'}');
    fclose(file_id);
    refused=false;
    try
        s2spinach(file_name);
    catch err
        refused=strcmp(err.identifier,bad_records{n,2})&&...
                contains(err.message,[file_name ':3:']);
    end
    result=test_true(result,bad_records{n,1},refused,...
                     'a diagnostic identifies both the invalid input and its source line.');

end

end


