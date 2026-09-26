% Tests phase-independent tensor-train compression error budgets. Syntax:
%
%                         result=test_tt_phase()
%
% Outputs:
%
%     result - regression result for signed and complex coefficients,
%              genuine rank reduction, zero boundaries, and spin products
%
% Finite tolerances bound the absolute Frobenius error. Zero tolerance
% comparisons allow numerical roundoff but require all resolved ranks.
%
% talos@spindynamics.org

function result=test_tt_phase()

% Initialise the regression record and reproducible complex fixtures
result=new_test_result('kernel/tt_phase','Tensor-train coefficient phases',...
                       'Compression must honour absolute error budgets independently of coefficient phase.');
saved_rng=rng; restore_rng=onCleanup(@()rng(saved_rng));
rng(240924); phases=[1 -1 1i -1i exp(0.37i) 2+3i];
max_error=0; reduced_cases=0;

% Exercise independent complex cores and ranks across two sweep lengths
for ncores=2:3
    cores=cell(ncores,3);
    for k=1:numel(cores)
        core=randn(2)+1i*randn(2);
        cores{k}=core/norm(core,'fro');
    end
    for nranks=1:3
        if nranks==1
            weights=1;
        elseif nranks==2
            weights=[1 1e-5];
        else
            weights=[1 0.15 1e-5];
        end

        % Keep phase in the coefficient of an already packed train
        P=pack(ttclass(weights,cores(:,1:nranks),1e-3*ones(1,nranks)/nranks));
        exact=pack(ttclass(weights,cores(:,1:nranks),zeros(1,nranks)));
        reference=full(P);
        for scale=[0.2 7]
            for phase=phases
                Q=(scale*phase)*P;
                rounded=shrink(Q);
                expected=(scale*phase)*reference;
                residual=norm(full(rounded)-expected,'fro');
                roundoff=1e-12*max(1,norm(expected,'fro'));
                max_error=max(max_error,residual);
                label=sprintf('%d cores, rank %d, scale %.1f, phase %.2f%+.2fi',...
                              ncores,nranks,scale,real(phase),imag(phase));
                result=test_close(result,label,full(rounded),expected,...
                                  Q.tolerance,1e-12,'absolute Frobenius truncation budget');

                % Require measurable truncation rather than an exact low-rank identity
                if nranks>1
                    assert(residual>100*eps*max(1,norm(expected,'fro')));
                    assert(all(rounded.ranks(2:end-1)<Q.ranks(2:end-1)));
                    reduced_cases=reduced_cases+1;
                end

                % Zero tolerance must retain the resolved rank and the dense value
                exact_round=shrink((scale*phase)*exact);
                assert(all(exact_round.ranks==exact.ranks));
                result=test_close(result,['zero tolerance: ' label],full(exact_round),...
                                  expected,roundoff,0,'zero tolerance permits only numerical roundoff');
            end
        end
    end
end

% Preserve the public compressor zero-coefficient escape at both tolerances
for tolerance=[0 1e-3]
    zero_train=shrink(ttclass(0,cores(:,1),tolerance));
    result=test_close(result,'zero coefficient',full(zero_train),zeros(2^ncores),...
                      0,0,'a zero coefficient must return an exact finite zero');
end

% Multiply a signed Hermitian spin coupling with genuinely complex local factors
spin_x=[0 1;1 0]/2; spin_y=[0 -1i;1i 0]/2; spin_z=[1 0;0 -1]/2;
left=(spin_x+spin_y)/sqrt(2); right=(spin_y+spin_z)/sqrt(2);
H=ttclass(-2*pi*150,{left;right},1e-8);
P=ttclass(1,{eye(2);eye(2)},0);
expected=-2*pi*150*kron(left,right);
assert(norm(expected-expected','fro')==0);
observed=full(H*P);
result=test_close(result,'signed Hermitian spin product',observed,expected,...
                  1e-8,0,'negative interaction coefficients must survive TT multiplication');
result=test_close(result,'spin product Hermiticity',observed,observed',...
                  1e-12,0,'compression must preserve the Hermitian product to roundoff');
fprintf('TT_PHASE max_error=%.15g reduced_cases=%d spin_error=%.15g\n',...
        max_error,reduced_cases,norm(observed-expected,'fro'));

end


