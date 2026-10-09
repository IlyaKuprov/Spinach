% Checks solvent-aware rigid-body surface-ellipsoid diffusion against
% sphere limits, solvent scaling, and coordinate invariants. Syntax:
%
%                         result=test_rotcorr()
%
% Outputs:
%
%    result - regression results and diagnostic messages
%
% talos@spindynamics.org

function result=test_rotcorr()

% Identify the physical contract
result=new_test_result('kernel/rotcorr','Surface ellipsoid diffusion',...
                      'Solvent models, sphere limits, and coordinate invariants.');

% Compare with the exact rank-2 Stokes-Einstein-Debye sphere limit
symbols={'H';'B';'C';'N';'O';'F';'Si';'P';'S';'Cl';'Se';'Br';'I'};
radii=[1.10;1.92;1.70;1.55;1.52;1.47;2.10;1.80;1.80;1.75;1.90;1.83;1.98];
visc=0.000538691125263246;
for n=1:numel(symbols)
    [tau,D,axes_len]=rotcorr(symbols(n),[0 0 0],'chloroform',298.15);
    reference=4*pi*visc*(radii(n)*1e-10)^3/(3*1.380649e-23*298.15);
    result=test_close(result,['sphere time ' symbols{n}],tau/reference,1,0,1e-6,...
                      'the sphere has tau=4*pi*eta*r^3/(3*k*T)');
    result=test_close(result,['sphere tensor ' symbols{n}],D*6*reference,...
                      eye(3),0,1e-4,'the sphere diffusion rates coincide');
    result=test_close(result,['sphere axes ' symbols{n}],axes_len/radii(n),...
                      ones(3,1),0,1e-4,'the bare sphere has the elemental radius');
end
[tau,~,axes_len]=rotcorr({'C'},[0 0 0],'water',298.15);
reference=4*pi*0.000889996773678783*(3.9e-10)^3/(3*1.380649e-23*298.15);
result=test_close(result,'hydrated sphere time',tau/reference,1,0,1e-6,...
                  'water adds a 2.2 Angstrom shell to the carbon radius');
result=test_close(result,'hydrated sphere axes',axes_len/3.9,ones(3,1),...
                  0,1e-4,'the probe is not added to the contact-surface radius');

% Exercise an anisotropic overlapping-sphere body
xyz=[-3 0 0;0 0 0;3 0 0]; atom_symbols={'C';'C';'C'};
[tau,D]=rotcorr(atom_symbols,xyz,'water',300);
result=test_true(result,'positive anisotropy',...
                 all(eig(D)>0)&&(max(eig(D))/min(eig(D))>1.1),...
                 'an elongated body has positive unequal diffusion rates');

% Check rigid translation without changing the surface quadrature
[moved,tensor]=rotcorr(atom_symbols,xyz+[17 -9 4],'water',300);
result=test_close(result,'translation time',moved/tau,1,0,1e-12,...
                  'translation cannot change rotational drag');
result=test_close(result,'translation tensor',tensor/norm(D),D/norm(D),...
                  0,1e-12,'translation cannot rotate the diffusion tensor');

% Check rigid rotation within the requested sampling tolerance
angle=0.71; R=[cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0;0 0 1];
[moved,tensor]=rotcorr(atom_symbols,xyz*R','water',300);
result=test_true(result,'rotation time',abs(moved/tau-1)<0.01,...
                 'rotated sampling changes tau by less than one percent');
result=test_true(result,'rotation tensor',...
                 max(abs(eig(tensor-R*D*R',tensor)))<0.01,...
                 'directional diffusion changes by less than one percent');

% Check that atom ordering does not change the hydrodynamic body
[moved,tensor]=rotcorr(atom_symbols,xyz([3 1 2],:),'water',300);
result=test_close(result,'permutation time',moved/tau,1,0,1e-12,...
                  'atom ordering cannot change rotational drag');
result=test_close(result,'permutation tensor',tensor/norm(D),D/norm(D),...
                  0,1e-12,'atom ordering cannot rotate the diffusion tensor');

% Check temperature scaling using the shipped viscosity model
solvents={'water','chloroform'};
viscosities=[0.00100156726460300 0.000693328713714806;...
             0.000566835797875055 0.000480862617762159];
for n=1:2
    cold=rotcorr(atom_symbols,xyz,solvents{n},293.15);
    warm=rotcorr(atom_symbols,xyz,solvents{n},310);
    ratio=viscosities(n,2)*293.15/(viscosities(n,1)*310);
    result=test_close(result,['temperature scaling ' solvents{n}],...
                      warm/cold,ratio,0,1e-12,...
                      'fixed-geometry tau scales as solvent viscosity over T');
end

% Check complete methane geometry in both solvent models
methane=[0 0 0;1 1 1;1 -1 -1;-1 1 -1;-1 -1 1]*1.09/sqrt(3);
for n=1:2
    [time,tensor,axes_len]=rotcorr({'C';'H';'H';'H';'H'},methane,...
                                  solvents{n},298.15);
    result=test_true(result,['methane ' solvents{n}],...
                     isfinite(time)&&time>0&&all(eig(tensor)>0)&&...
                     all(axes_len>0),...
                     'a complete tetrahedral molecule has positive finite outputs');
end

% Reject unsupported chemistry and malformed physical inputs
bad_inputs={{'Xx'},[0 0 0],'water',300;...
            {'C'},[0 0 0],"water",300;...
            {'C'},[0 0 0],'ethanol',300;...
            {'C'},[0 0 0],'water',0;...
            {'C'},[0 0 0],'chloroform',1000;...
            {'C'},[0 NaN 0],'water',300;...
            {'C','C'},[0 0 0;1 0 0],'water',300;...
            {'C';'C'},[0 0 0;0 0 0],'water',300};
for n=1:size(bad_inputs,1)
    rejected=false;
    try
        rotcorr(bad_inputs{n,:});
    catch
        rejected=true;
    end
    result=test_true(result,['invalid input ' num2str(n)],rejected,...
                     'unsupported or malformed inputs must raise an error');
end

end


