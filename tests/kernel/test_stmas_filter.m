% Checks the STMAS transition support projector for spin 3/2 and
% spin 5/2, including spectators on both sides of the selected spin.
% Syntax:
%
%                     result=test_stmas_filter()
%
% Outputs:
%
%     result - regression test result with explanatory messages
%
% talos@spindynamics.org

function result=test_stmas_filter()

% State the target of the test
result=new_test_result('kernel/transition_filter',...
                       'STMAS transition support projector',...
                       'selected Zeeman coherences retain their amplitudes and spectator states.');

% Use a quadrupolar spin between unlike spectators in product order
sys.magnet=0;
sys.output='hush';
sys.parallel={'processes',1};
inter.zeeman.scalar={0,0,0};
bas.formalism='zeeman-liouv';
bas.approximation={'none'};
isotopes={'87Rb','27Al'};
for n=1:numel(isotopes)

    % Build the system and identify central and satellite level pairs
    sys.isotopes={'1H',isotopes{n},'2H'};
    spin_system=basis(create(sys,inter),bas);
    mult=spin_system.comp.mults(2);
    ct_pair=[mult/2 mult/2+1];
    rows=(1:mult-1)';
    sat_pairs=[rows rows+1];
    sat_pairs(ct_pair(1),:)=[];
    label=[isotopes{n} ' '];

    % Reconstruct the satellite support after subtracting the central ladder element
    ladder=diag(sqrt(rows.*(mult-rows)),1);
    ct=centrans(mult,'+');
    [sat_rows,sat_cols]=find(ladder-(mult/2)*ct);
    result=test_close(result,[label 'satellite branches'],[sat_rows sat_cols],sat_pairs,0,0,...
                      'both satellite branches are in the positive coherence order and exclude the CT');

    % Compare positive satellite selection against explicit Hilbert-space projectors
    dim=6*mult;
    rho=reshape((1:dim^2).'+1i*(dim^2:-1:1).',dim,dim);
    retained=zeros(dim);
    for k=1:size(sat_pairs,1)

        % Form each term P_a rho P_b with spectator identities
        row_proj=sparse(sat_pairs(k,1),sat_pairs(k,1),1,mult,mult);
        col_proj=sparse(sat_pairs(k,2),sat_pairs(k,2),1,mult,mult);
        row_proj=kron(speye(2),kron(row_proj,speye(3)));
        col_proj=kron(speye(2),kron(col_proj,speye(3)));
        retained=retained+row_proj*rho*col_proj;

    end
    P=stmas_filter(spin_system,2,sat_pairs);
    result=test_close(result,[label 'selected entries'],P*rho(:),retained(:),0,0,...
                      'complex amplitudes are unchanged, including off-diagonal spectator entries');
    result=test_true(result,[label 'sparse idempotence'],issparse(P)&&isequal(P*P,P),...
                     'a second application neither scales nor removes retained entries');
    result=test_true(result,[label 'duplicate integer pairs'],...
                     isequal(P,stmas_filter(spin_system,uint8(2),uint8([sat_pairs; sat_pairs]))),...
                     'integer level indices are accepted and duplicate pairs form a union');

    % Check the CT normalisation and directional isolation in the native representation
    ct_full=kron(speye(2),kron(ct,speye(3)));
    ct_state=coil_state(spin_system,'CT+',2);
    result=test_close(result,[label 'unit CT amplitude'],ct_state,ct_full(:),1e-13,0,...
                      'Spinach CT+ is a unit matrix element, not the central ladder coefficient');
    P_ct=stmas_filter(spin_system,2,ct_pair);
    result=test_close(result,[label 'CT isolation'],P_ct*ct_state,ct_full(:),1e-13,0,...
                      'the CT projector preserves the central coherence');
    result=test_true(result,[label 'rejected coherences'],nnz(P_ct*P)==0&&...
                     nnz(P_ct*coil_state(spin_system,'CT-',2))==0,...
                     'central and satellite coherences are disjoint, and reverse CT is rejected');

end

% Refuse level indices outside the selected spin
try
    stmas_filter(spin_system,2,[1 mult+1]);
    refused=false;
catch err
    refused=strcmp(err.identifier,'Spinach:stmas_filter:levels');
end
result=test_true(result,'invalid level index',refused,...
                 'a level outside the selected spin is reported rather than indexed');

% Refuse a different coordinate representation
bas.formalism='zeeman-hilb';
spin_system=basis(spin_system,bas);
try
    stmas_filter(spin_system,2,ct_pair);
    refused=false;
catch err
    refused=strcmp(err.identifier,'Spinach:stmas_filter:basis');
end
result=test_true(result,'invalid formalism',refused,...
                 'column-stretched masks do not act directly on Hilbert matrices');

% An empty selection is the zero projector
bas.formalism='zeeman-liouv';
spin_system=basis(spin_system,bas);
P=stmas_filter(spin_system,2,zeros(0,2));
result=test_true(result,'empty selection',nnz(P)==0&&isequal(size(P),[dim^2 dim^2]),...
                 'an empty union keeps no elements but preserves the Liouville dimension');

end


