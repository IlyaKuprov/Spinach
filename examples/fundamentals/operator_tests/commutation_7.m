% Expansion relations for operator basis transforms.
%
% ilya.kuprov@weizmann.ac.il

function commutation_7()

% Accuracy threshold
tol=1e-10;

% Test IST expansion of Zeeman level projectors
mult=5; T=irr_sph_ten(mult);
for lvl_num=1:mult

    % Build known Zeeman projector
    P=zeros(mult,mult);
    P(mult-lvl_num+1,mult-lvl_num+1)=1;

    % Obtain IST expansion.
    [states,coeffs]=enlev2ist(mult,lvl_num);

    % Reconstruct operator from IST terms
    P_rec=zeros(mult,mult,'like',1i);
    for n=1:numel(states)
        P_rec=P_rec+coeffs(n)*T{states(n)+1};
    end

    % Report IST expansion failures
    if norm(P-P_rec,'fro')>tol
        error('Zeeman projector IST expansion test FAILED.');
    end

end
disp('Zeeman projector IST expansion test PASSED.');

% Test BM expansion of bosonic level projectors
nlevels=6; B=boson_mono(nlevels);
for lvl_num=1:nlevels

    % Build known bosonic projector
    P=zeros(nlevels,nlevels);
    P(lvl_num,lvl_num)=1;

    % Obtain BM expansion.
    [states,coeffs]=enlev2bm(nlevels,lvl_num);

    % Reconstruct operator from BM terms
    P_rec=zeros(nlevels,nlevels,'like',1i);
    for n=1:numel(states)
        P_rec=P_rec+coeffs(n)*B{states(n)+1};
    end

    % Report BM expansion failures
    if norm(P-P_rec,'fro')>tol
        error('Bosonic projector BM expansion test FAILED.');
    end

end
disp('Bosonic projector BM expansion test PASSED.');

% Test bosonic product to IST conversion
W=weyl(nlevels); T_bos=irr_sph_ten(nlevels);
prod_list={'CA','ACCA','CCAAA'};
for n=1:numel(prod_list)

    % Build operator product explicitly
    P=speye(nlevels,nlevels);
    for k=1:numel(prod_list{n})
        if prod_list{n}(k)=='C'
            P=P*W.c;
        else
            P=P*W.a;
        end
    end

    % Obtain IST expansion
    [states,coeffs]=bos2ist(prod_list{n},nlevels);

    % Reconstruct operator from IST terms
    P_rec=zeros(nlevels,nlevels,'like',1i);
    for k=1:numel(states)
        P_rec=P_rec+coeffs(k)*T_bos{states(k)+1};
    end

    % Report bosonic IST conversion failures
    if norm(P-P_rec,'fro')>tol
        error('Bosonic operator IST conversion test FAILED.');
    end

end
disp('Bosonic operator IST conversion test PASSED.');

% Test matrix expansion into single-transition basis
dim=7; A=randn(dim)+1i*randn(dim);
E=sin_tran(dim); A_rec=zeros(dim,dim,'like',1i);
for n=1:numel(E)

    % Add transition basis term
    A_rec=A_rec+hdot(E{n},A)*E{n};

end

% Report single-transition expansion failures
if norm(A-A_rec,'fro')>tol
    error('Single transition basis expansion test FAILED.');
end
disp('Single transition basis expansion test PASSED.');

% Test BM basis expansion of a random bosonic operator
A=randn(nlevels)+1i*randn(nlevels);
[states,coeffs]=oper2bm(A);
A_rec=zeros(nlevels,nlevels,'like',1i);
for n=1:numel(states)

    % Add BM basis term
    A_rec=A_rec+coeffs(n)*B{states(n)+1};

end

% Report BM expansion failures
if norm(A-A_rec,'fro')/norm(A,'fro')>1e-8
    error('Bosonic monomial expansion test FAILED.');
end
disp('Bosonic monomial expansion test PASSED.');

end

