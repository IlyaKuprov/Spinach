% Stream the secular spin-1/2 Hahn pair product. Syntax:
%
%  [echo,omitted]=cce_pair_product(xyz,axis,hf,dipolar,exchange,times,atol,bytes)
%
% Parameters:
%
%   xyz     - N-by-3 nuclear coordinates in a common length unit
%   axis    - 3-by-1 unit field direction in the coordinate frame
%   hf      - N-by-1 longitudinal hyperfine coefficients in rad/s
%   dipolar - nuclear point-dipole prefactor in rad/s times length cubed
%   exchange  - sparse N-by-N lower triangle of isotropic bath exchange in rad/s
%   times   - nonnegative column of total echo times in seconds
%   atol    - absolute echo-loss budget for screening weak pairs
%   bytes   - positive working-memory budget in bytes
%
% Outputs:
%
%   echo    - column of normalised Hahn echo amplitudes
%   omitted - uniform bound on screened echo loss, excluding roundoff
%
% The two branches in the flip-flop subspace are
% H_±=[-(Dzz-2*J)*sigma_x ± delta_A*sigma_z]/4, up to identity.
%
% ilya.kuprov@weizmann.ac.il

function [echo,omitted]=cce_pair_product(xyz,axis,hf,dipolar,exchange,times,atol,bytes)

% Check consistency
grumble(xyz,axis,hf,dipolar,exchange,times,atol,bytes);

% Budget row work arrays and simultaneous pair-by-time temporaries
n_spins=size(xyz,1);
block=min(max(n_spins,1),floor(bytes/(32*(4+numel(times)))));
if block<1
    error('cce_pair_product:memory','The memory budget cannot hold one pair trace.');
end
pair_tol=atol/max(1,n_spins*(n_spins-1)/2);
angles=times.'/8;
max_angle=max(angles);
log_echo=zeros(1,numel(times));
omitted=0;
has_exchange=nnz(exchange)>0;
for n=1:(n_spins-1)
    for k=(n+1):block:n_spins

        % Compute only this row block of pair interactions
        partners=k:min(k+block-1,n_spins);
        delta=xyz(partners,:)-xyz(n,:);
        dist=vecnorm(delta,2,2);
        pair=dipolar*(1-3*((delta*axis)./dist).^2)./dist.^3;
        if has_exchange
            pair=pair-2*full(exchange(partners,n));
        end
        detuning=hf(n)-hf(partners);

        % Zero splittings refocus exactly; contrast is bounded by one
        freq=hypot(pair,detuning);
        if any(~isfinite(freq))
            error('cce_pair_product:splitting','Pair splittings must be finite.');
        end
        active=freq>0;
        freq=freq(active);
        contrast=min(1,(2*(pair(active)./freq).*(detuning(active)./freq)).^2);

        % Bound weak-pair losses by sin(x)^4<=min(1,x^4)
        bound=contrast.*min(1,(max_angle*freq).^4);
        keep=bound>pair_tol;
        omitted=omitted+sum(bound(~keep));
        if any(keep)
            losses=contrast(keep).*sin(freq(keep).*angles).^4;
            log_echo=log_echo+sum(log1p(-losses),1);
        end

    end
end
echo=exp(log_echo).';

end

% Check the numerical kernel inputs and dimensions
function grumble(xyz,axis,hf,dipolar,exchange,times,atol,bytes)
n_spins=size(xyz,1);
if ~isnumeric(xyz)||~isreal(xyz)||~ismatrix(xyz)||size(xyz,2)~=3||...
   any(~isfinite(xyz),'all')||size(unique(xyz,'rows'),1)~=n_spins
    error('xyz must be an N-by-3 real array of finite distinct sites.');
end
if ~isnumeric(axis)||~isreal(axis)||~isequal(size(axis),[3 1])||...
   ~(abs(norm(axis)-1)<=4*eps)
    error('axis must be a three-element real unit column vector.');
end
if ~isnumeric(hf)||~isreal(hf)||~isequal(size(hf),[n_spins 1])||any(~isfinite(hf))
    error('hf must be a finite real N-by-1 array.');
end
if ~isnumeric(dipolar)||~isreal(dipolar)||~isscalar(dipolar)||~isfinite(dipolar)
    error('dipolar must be a finite real scalar.');
end
if ~issparse(exchange)||~isreal(exchange)||~isequal(size(exchange),[n_spins n_spins])
    error('exchange must be a real N-by-N sparse matrix.');
end
if ~isnumeric(times)||~isreal(times)||~iscolumn(times)||isempty(times)||...
   any(~isfinite(times)|times<0)
    error('times must be a nonempty column of finite nonnegative values.');
end
if ~isnumeric(atol)||~isreal(atol)||~isscalar(atol)||~isfinite(atol)||atol<0
    error('atol must be a finite nonnegative real scalar.');
end
if ~isnumeric(bytes)||~isreal(bytes)||~isscalar(bytes)||~(bytes>0)
    error('bytes must be a positive real scalar.');
end
end


