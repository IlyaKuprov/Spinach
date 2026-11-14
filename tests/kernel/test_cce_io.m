% Versioned complex result round trips and explicit physical input units.
%
%                         result=test_cce_io()
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_io()

% Make isolated temporary files and standard harness output
result=new_test_result('kernel/cce_io','CCE input and output',...
                       'Complex shape/type round trips, scalar/tensor units and cubes.');
folder=tempname; mkdir(folder); cleanup=onCleanup(@()rmdir(folder,'s'));
original=struct('value',[1+2i NaN;3-4i Inf],'validity',logical([1 0;1 0]),...
                'diagnostics',struct('bounds',[eps 2*eps],'reason',{{'valid','zero'}}),...
                'empty',zeros(0,2),'tensor',reshape(1:24,2,3,4),...
                'id',intmax('uint64'),'powers',[intmin('int64') intmax('int64')]);

% Test each versioned representation, including its complete diagnostics
for format={'mat','json','csv'}
    file=fullfile(folder,['result.' format{1}]); cce_export(original,file,format{1});
    restored=cce_import(file,format{1},struct());
    result=test_true(result,['Exact complex round trip: ' format{1}],isequaln(original,restored),['Exact complex round trip: ' format{1}]);
end

% Verify declared coordinate units and active Cartesian rotations
file=fullfile(folder,'sites.xyz'); put(file,sprintf('2\nexplicit isotope labels\n1H 0.1 0.2 0.3\n13C 0.4 0.5 0.6\n'));
u=struct('length','nm','frame','molecular','rotation',[0 -1 0;1 0 0;0 0 1]);
data=cce_import(file,'xyz',u);
result=test_close(result,'XYZ units and active rotation',data.coordinates_A,[-2 1 3;-5 4 6],1e-14,1e-14,'XYZ units and active rotation');

% Compare scalar and full tensor couplings under independent unit conversions
file=fullfile(folder,'scalar.csv'); put(file,sprintf('1,2,2\n'));
u=struct('tensor','MHz','frame','lab'); data=cce_import(file,'tensor',u);
result=test_close(result,'Scalar frequency is isotropic tensor',data.pairs_Hz.tensor,2e6*eye(3),0,0,'Scalar frequency is isotropic tensor');
put(file,sprintf('1,2,1,2,3,4,5,6,7,8,9\n')); u.tensor='rad/s';
data=cce_import(file,'tensor',u);
result=test_close(result,'Row-major tensor and angular frequency',data.pairs_Hz.tensor,[1 2 3;4 5 6;7 8 9]/(2*pi),1e-15,1e-15,'Row-major tensor and angular frequency');
u.tensor='G'; u.gyro_rad_s_T=2*pi*1e6; put(file,sprintf('1,2,2\n'));
data=cce_import(file,'tensor',u);
result=test_close(result,'Field units require an explicit gyro',data.pairs_Hz.tensor,200*eye(3),1e-12,1e-14,'Field units require an explicit gyro');
u=rmfield(u,'gyro_rad_s_T'); rejected=false;
try, cce_import(file,'tensor',u); catch, rejected=true; end
result=test_true(result,'Unspecified field gyromagnetic units are refused',rejected,'Unspecified field gyromagnetic units are refused');

% Check cube geometry and z-fastest voxel ordering against an explicit array
file=fullfile(folder,'density.cube');
put(file,sprintf('density\nminimal\n1 0 0 0\n2 1 0 0\n1 0 1 0\n2 0 0 1\n1 0 0 0 0\n1 2 3 4\n'));
u=struct('length','bohr','frame','lab','density','electron/bohr^3');
data=cce_import(file,'cube',u);
result=test_close(result,'Cube z-fastest order',squeeze(data.density),[1 2;3 4],0,0,'Cube z-fastest order');
result=test_close(result,'Cube explicit bohr conversion',data.voxel_vectors_A,0.529177210544*eye(3),0,0,'Cube explicit bohr conversion');
disp('CCE IO: MAT/JSON/CSV round trips exact; tensor and voxel conversions passed.');

end

function put(file,content)

% Write a bounded test input and close it immediately
fid=fopen(file,'w'); cleanup=onCleanup(@()fclose(fid)); fwrite(fid,content,'char');

end
