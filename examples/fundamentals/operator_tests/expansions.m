% Product action tests for irreducible spherical tensor
% operators and orthogonalised bosonic monomials.
%
% ilya.kuprov@weizmann.ac.il

function expansions()

% Random level counts
spin_mult=1+randi(9);
bos_nlevels=2+randi(8);

% Monomials and ISTs
T=irr_sph_ten(spin_mult);
B=boson_ortho(bos_nlevels);

% Multiplication tables
[ist_PTL,ist_PTR]=ist_product_table(spin_mult);
[bos_PTL,bos_PTR]=bos_product_table(bos_nlevels);

% IST test loop
for n=1:numel(T)
    for m=1:numel(T)

        % Left side product action
        LSPA_table=zeros(spin_mult,spin_mult);
        for k=1:numel(T)
            LSPA_table=LSPA_table+ist_PTL(n,m,k)*T{k}/norm(T{k},'fro');
        end
        LSPA_known=T{n}*T{m}/norm(T{m},'fro');

        % Right side product action
        RSPA_table=zeros(spin_mult,spin_mult);
        for k=1:numel(T)
            RSPA_table=RSPA_table+ist_PTR(n,m,k)*T{k}/norm(T{k},'fro');
        end
        RSPA_known=T{m}*T{n}/norm(T{m},'fro');

        % Accuracy tests
        norm_diff_a=norm(LSPA_table-LSPA_known,'fro');
        norm_diff_b=norm(RSPA_table-RSPA_known,'fro');
        norm_scale=norm(T{m},'fro'); tol=sqrt(eps);
        if (norm_diff_a/norm_scale>tol)||...
           (norm_diff_b/norm_scale>tol)
            error('IST product action table test FAILED.');
        end

    end
end
disp('IST product action table test PASSED.');

% BM test loop
for n=1:numel(B)
    for m=1:numel(B)

        % Left side product action
        LSPA_table=zeros(bos_nlevels,bos_nlevels);
        for k=1:numel(B)
            LSPA_table=LSPA_table+bos_PTL(n,m,k)*B{k}/norm(B{k},'fro');
        end
        LSPA_known=B{n}*B{m}/norm(B{m},'fro');

        % Right side product action
        RSPA_table=zeros(bos_nlevels,bos_nlevels);
        for k=1:numel(B)
            RSPA_table=RSPA_table+bos_PTR(n,m,k)*B{k}/norm(B{k},'fro');
        end
        RSPA_known=B{m}*B{n}/norm(B{m},'fro');

        % Accuracy tests
        norm_diff_a=norm(LSPA_table-LSPA_known,'fro');
        norm_diff_b=norm(RSPA_table-RSPA_known,'fro');
        norm_scale=norm(B{m},'fro'); tol=sqrt(eps);
        if (norm_diff_a/norm_scale>tol)||...
           (norm_diff_b/norm_scale>tol)
            error('BM product action table test FAILED.');
        end

    end
end
disp('BM product action table test PASSED.');

end

