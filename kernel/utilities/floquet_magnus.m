% Fourier-space Floquet-Magnus expansion of a periodic matrix generator.
% Computes the effective generator and periodic kick without time integrals,
% following Appendix A of https://doi.org/10.5194/mr-2026-10. Syntax:
%
%          [F,K]=floquet_magnus(H,omega,order,t0,theory)
%
% Parameters:
%
%     H     - row cell array {H_-N,...,H_0,...,H_N} of equally sized
%             square double matrices, full or sparse, in rad/s;
%             H(t)=sum_n H_n*exp(1i*n*omega*t)
%
%     omega - finite nonzero real double scalar, angular frequency in rad/s
%
%     order - real double scalar, integer inverse-frequency power 0 to 3;
%             order 0 is H_0, order 3 is fourth order in the paper
%
%     t0    - finite real double scalar, initial time in seconds
%
%     theory - 'van_vleck' for zero-mean kick coefficients, or 'aht'
%              for a kick that vanishes at t0 order by order
%
% Outputs:
%
%     F     - effective generator through omega^(-order), rad/s
%
%     K     - row cell array of dimensionless kick Fourier coefficients
%             with harmonic indices -(order+1)*N:(order+1)*N, including
%             all powers through omega^(-order-1). Evaluate this series
%             as K(t)=sum_n K_n*exp(1i*n*omega*t). The approximation is
%             U(t,t0)=expm(-1i*K(t))*expm(-1i*F*(t-t0))*expm(1i*K(t0)).
%
% Notes: this is a high-frequency asymptotic expansion, not a convergence
% guarantee. Retaining one extra kick order supports intra-cycle sampling.
% AHT and Van Vleck truncations agree asymptotically, not exactly. Hermitian
% Hamiltonians have H_-n=H_n'; then F and K(t) are Hermitian up to roundoff.
% General matrix inputs are allowed, but do not imply unitary evolution,
% positivity preservation, or a relaxation-theory approximation.
%
% talos@spindynamics.org

function [F,K]=floquet_magnus(H,omega,order,t0,theory)

% Check the caller-controlled inputs
grumble(H,omega,order,t0,theory);

% Retain the full Fourier support generated at the requested order
in_rank=(numel(H)-1)/2; out_rank=(order+1)*in_rank;
centre=out_rank+1; ncoeff=2*out_rank+1; dim=size(H{1},1);
zero=repmat({sparse(dim,dim)},1,ncoeff);
harmonics=zero; harmonics(centre-in_rank:centre+in_rank)=H;
kick_terms=repmat(zero,order+1,1); eff_terms=cell(1,order+1);

% Solve the Fourier recurrence in increasing inverse-frequency power
for power=0:order

    % The zeroth coefficient is the supplied time-dependent generator
    if power==0
        rhs=harmonics;
    else
        rhs=zero;

        % Single commutators carry the first Bernoulli coefficient
        term=fourier_comm(kick_terms(power,:),harmonics);
        for n=1:ncoeff
            rhs{n}=rhs{n}+(1i/2)*term{n};
        end
        for a=1:power
            constant=zero; constant{centre}=eff_terms{power-a+1};
            term=fourier_comm(kick_terms(a,:),constant);
            for n=1:ncoeff
                rhs{n}=rhs{n}+(1i/2)*term{n};
            end
        end

        % Double commutators carry B_2/2! and the opposite F sign
        for a=1:power-1
            for b=1:power-a
                constant=zero; constant{centre}=-eff_terms{power-a-b+1};
                if a+b==power
                    for n=1:ncoeff
                        constant{n}=constant{n}+harmonics{n};
                    end
                end
                term=fourier_comm(kick_terms(b,:),constant);
                term=fourier_comm(kick_terms(a,:),term);
                for n=1:ncoeff
                    rhs{n}=rhs{n}-term{n}/12;
                end
            end
        end

    end

    % The zero Fourier component supplies the next effective generator
    eff_terms{power+1}=rhs{centre};

    % Integrate only nonzero harmonics to obtain the next periodic kick
    for n=[1:centre-1 centre+1:ncoeff]
        kick_terms{power+1,n}=rhs{n}/(1i*(n-centre)*omega);
    end

    % Fix the integration constant in the requested gauge
    if strcmp(theory,'aht')
        for n=[1:centre-1 centre+1:ncoeff]
            kick_terms{power+1,centre}=kick_terms{power+1,centre}-...
                                       kick_terms{power+1,n}*exp(1i*(n-centre)*omega*t0);
        end
    end

end

% Sum the explicitly computed orders without altering their symmetry
F=sparse(dim,dim); K=zero;
for power=1:order+1
    F=F+eff_terms{power};
    for n=1:ncoeff
        K{n}=K{n}+kick_terms{power,n};
    end
end

end

% Fourier convolution of a matrix commutator without cyclic wraparound
function C=fourier_comm(A,B)

% The allocated support contains every harmonic needed at this order
ncoeff=numel(A); centre=(ncoeff+1)/2;
C=repmat({sparse(size(A{1},1),size(A{1},2))},1,ncoeff);

% Accumulate nonzero matrix products at the sum of harmonic indices
for a=1:ncoeff
    if nnz(A{a})==0, continue; end
    for b=max(1,centre+1-a):min(ncoeff,ncoeff+centre-a)
        if nnz(B{b})==0, continue; end
        n=a+b-centre;
        C{n}=C{n}+A{a}*B{b}-B{b}*A{a};
    end
end

end

% Validate the Fourier input and expansion controls
function grumble(H,omega,order,t0,theory)
if ~iscell(H)||~isrow(H)||isempty(H)||mod(numel(H),2)~=1
    error('H must be a nonempty odd-length row cell array.');
end
if ~isa(H{1},'double')||~ismatrix(H{1})||isempty(H{1})||...
   size(H{1},1)~=size(H{1},2)
    error('Fourier coefficients must be nonempty square double matrices.');
end
for n=1:numel(H)
    if ~isa(H{n},'double')||~ismatrix(H{n})||...
       ~isequal(size(H{n}),size(H{1}))||any(~isfinite(H{n}),'all')
        error('Fourier coefficients must be equally sized finite double matrices.');
    end
end
if ~isa(omega,'double')||~isreal(omega)||~isscalar(omega)||...
   ~isfinite(omega)||omega==0
    error('omega must be a finite nonzero real double scalar.');
end
if ~isa(order,'double')||~isreal(order)||~isscalar(order)||...
   ~isfinite(order)||order~=fix(order)||order<0||order>3
    error('order must be a double scalar integer between 0 and 3.');
end
if ~isa(t0,'double')||~isreal(t0)||~isscalar(t0)||~isfinite(t0)
    error('t0 must be a finite real double scalar.');
end
if ~ischar(theory)||~isrow(theory)||...
   ~ismember(theory,{'aht','van_vleck'})
    error('theory must be ''aht'' or ''van_vleck''.');
end
end


