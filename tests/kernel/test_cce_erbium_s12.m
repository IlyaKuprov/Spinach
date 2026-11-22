% Test: kernel/cce_erbium_s12 | Er:CaWO4 secular pair reduction
% Compare both analytical Er:CaWO4 CCE-2 routes to dense, full pair
% Hahn propagation. Returns the standard Spinach test result record.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_erbium_s12()

% Three tungsten sites in the region of resolvable flip-flop decay
result=new_test_result('kernel/cce_erbium_s12','Er:CaWO4 secular pair reduction',...
                       'Compare the analytical pair product with explicit Hahn propagation.');
xyz=1e-9*[4.456465 0.000000 -0.2843425;...
          4.194320 0.262145  0.2843425;...
          4.456465 0.262145  0.0000000];
field_axis=[cosd(46.5);sind(46.5);0];
gamma_w=spin('183W');
hyperfine=1e-7*8.38*9.2740100657e-24*gamma_w;
dipolar=1e-7*(6.62607015e-34/(2*pi))*gamma_w^2;
times=linspace(0,0.060,61).';

% Construct full four-dimensional secular pair Hamiltonians
sx=[0 1;1 0]/2;
sy=[0 -1i;1i 0]/2;
sz=diag([1 -1])/2;
pair_op=kron(sz,sz)-(kron(sx,sx)+kron(sy,sy))/2;
zeeman=-gamma_w*0.067*(kron(sz,eye(2))+kron(eye(2),sz));
reference=ones(size(times));
pair_ref=ones(size(times));
for j=1:size(xyz,1)-1
    for k=j+1:size(xyz,1)

        % Compute dipolar tensor projections independently of the analytical function
        r=xyz(k,:)-xyz(j,:);
        u=r.'/norm(r);
        D=dipolar*(eye(3)-3*(u*u.'))/norm(r)^3;
        r=xyz(j,:);
        u=r.'/norm(r);
        A=hyperfine*(eye(3)-3*(u*u.'))/norm(r)^3;
        r=xyz(k,:);
        u=r.'/norm(r);
        B=hyperfine*(eye(3)-3*(u*u.'))/norm(r)^3;
        coupling=(field_axis.'*A*field_axis)*kron(sz,eye(2))+...
                 (field_axis.'*B*field_axis)*kron(eye(2),sz);
        H=zeeman+(field_axis.'*D*field_axis)*pair_op;
        for n=1:numel(times)

            % Compute the trace of the full Hahn commutator including spectator levels
            ua=expm(-1i*(H+coupling/2)*times(n)/2);
            ub=expm(-1i*(H-coupling/2)*times(n)/2);
            factor=trace((ua*ub)'*(ub*ua))/4;
            reference(n)=reference(n)*factor;
            if j==1&&k==2
                pair_ref(n)=factor;
            end

        end

    end
end

% Compare the analytical product to the explicit propagators
answer=cce_secular_pairs(xyz,field_axis,hyperfine,dipolar,times,0);
result=test_true(result,'Resolved decay',min(real(reference))<0.999,...
                 'The finite bath must distinguish the echo from unity.');
result=test_close(result,'Secular pair product',answer,reference,1e-8,1e-8,...
                  'Analytical flip-flop blocks must match dense propagation.');

% A singleton and a commuting nuclear bath both have perfect Hahn refocusing
for k=[1 size(xyz,1)]

    % Vanishing detuning exercises rows with no retained factors
    answer=cce_secular_pairs(xyz(1:k,:),field_axis,0,dipolar,times,0);
    result=test_true(result,sprintf('Commuting case (%d spins)',k),all(answer==1),...
                     'Without hyperfine detuning the Hahn echo must equal unity.');

end

% Create the same system from coordinates with the dipolar network switched off
sys.magnet=0.067;
sys.output='hush';
sys.parallel={'processes',1};
sys.tols.prox_cutoff=0;
parameters.observed=1;
parameters.transition=struct('spins',1,'ms',[1/2 -1/2]);
parameters.orientation=deg2rad([0 90 133.5]);
parameters.times=times.';
parameters.cce.order=2;
parameters.cce.approximation='projected';
parameters.cce.backend='analytical_hahn';
parameters.cce.numerics=struct('atol',0,'max_memory_bytes',32*(4+numel(times)));

% Exercise both a single pair and the complete three-nucleus pair family
expected=[pair_ref reference];
for n=2:3
    sys.isotopes=[{'E'} repmat({'183W'},1,n)];
    inter.zeeman.scalar=cell(size(sys.isotopes));
    inter.zeeman.scalar{1}=8.38;
    inter.coordinates=num2cell([zeros(1,3); 1e10*xyz(1:n,:)],2);
    spin_system=create(sys,inter);
    [answer,details]=clexpn(spin_system,@hahn_echo,parameters,'labframe');
    result=test_close(result,sprintf('Analytical Hahn with %d nuclei',n),...
                      answer,expected(:,n-1),1e-8,1e-8,...
                      'Streamed pair factors must match full pair propagation.');
    result=test_true(result,'Exact singleton and observed-only echoes',...
                     all(details.reduction.by_order(:,:,1:2)==1,'all'),...
                     'The reference and every singleton refocus exactly.');
end

% Opposite explicit hyperfine tensors cancel the native electron dipole
mixed=inter;
mixed.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes));
for k=1:size(xyz,1)
    r=xyz(k,:); u=r.'/norm(r);
    mixed.coupling.matrix{1,k+1}=hyperfine*(eye(3)-3*(u*u.'))/(2*pi*norm(r)^3);
end
mixed_system=create(sys,mixed);
cancelled=clexpn(mixed_system,@hahn_echo,parameters,'labframe');
result=test_close(result,'Mixed hyperfine cancellation',cancelled,ones(size(times)),1e-12,0,...
                  'Measured and point-dipole tensors must add with the native electron sign.');

% Unequal nuclear shifts do not refocus as a common frequency
mixed.zeeman.scalar{2}=10;
mixed_system=create(sys,mixed);
rejected=false;
try
    clexpn(mixed_system,@hahn_echo,parameters,'labframe');
catch exception
    rejected=contains(exception.message,'common nuclear Zeeman');
end
result=test_true(result,'Unequal nuclear shifts',rejected,'A bath detuning must not be silently discarded.');

% The reported pruning bound must contain the measured change from no screening
exact=answer;
parameters.cce.numerics.atol=3;
[answer,details]=clexpn(spin_system,@hahn_echo,parameters,'labframe');
bound=details.reduction.abs_error_bound;
result=test_true(result,'Reported screening bound',...
                 any(bound>0)&&all(abs(answer-exact)<=bound+eps),...
                 'The uniform pair-loss bound must cover the screened result.');
parameters.cce.numerics.atol=0;

% With no couplings every flip-flop splitting vanishes
inter.coordinates{1}=[];
spin_system=create(sys,inter);
answer=clexpn(spin_system,@hahn_echo,parameters,'labframe');
result=test_true(result,'Zero splitting',all(answer==1),...
                 'An exactly degenerate pair contributes unity without a division.');

% Bath species and couplings outside the spin-1/2 exchange model must be rejected
deuterons=sys;
deuterons.isotopes(2:end)={'2H'};
skewed=inter;
skewed.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes));
skewed.coupling.matrix{2,3}=[0 1 0; -1 0 0; 0 0 0];
cases={create(deuterons,inter),create(sys,skewed)};
names={'Spin-1 bath rejected','Antisymmetric bath tensor rejected'};
for k=1:2

    % Require the evaluator's own error rather than a silent unit echo
    rejected=false;
    try
        clexpn(cases{k},@hahn_echo,parameters,'labframe');
    catch exception
        rejected=strcmp(exception.identifier,'clexpn:analyticalHahn');
    end
    result=test_true(result,names{k},rejected,...
                     'The pair formula covers spin-1/2 nuclei with isotropic exchange only.');

end

% Reject requests whose physics the analytical evaluator does not implement
parameters.cce.graph=struct('mode','distance','radius',10);
rejected=false;
try
    clexpn(spin_system,@hahn_echo,parameters,'labframe');
catch exception
    rejected=strcmp(exception.identifier,'clexpn:analyticalHahn');
end
result=test_true(result,'Unsupported graph rejected',rejected,...
                 'Complete-family evaluation cannot silently change graph policy.');

end


