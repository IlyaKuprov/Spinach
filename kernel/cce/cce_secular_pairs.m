% Secular CCE-2 Hahn echo for an unpolarised bath of identical spin-1/2
% nuclei with longitudinal dipolar hyperfine couplings. Syntax:
%
%  echo=cce_secular_pairs(xyz,field_axis,hyperfine,dipolar,times,abs_tol)
%
% Parameters:
%
%     xyz        - N-by-3 array of centre-relative nuclear sites, in metres
%
%     field_axis - 3-by-1 magnetic field unit vector
%
%     hyperfine  - point-dipole hyperfine prefactor, in rad s^-1 m^3
%
%     dipolar    - nuclear point-dipole prefactor, in rad s^-1 m^3
%
%     times      - column of total Hahn echo times, in seconds
%
%     abs_tol    - upper bound on echo loss omitted by pair screening
%
% Outputs:
%
%     echo       - column of dimensionless normalised echo amplitudes
%
% ilya.kuprov@weizmann.ac.il

function echo=cce_secular_pairs(xyz,field_axis,hyperfine,dipolar,times,abs_tol)

% Check consistency
grumble(xyz,field_axis,hyperfine,dipolar,times,abs_tol);

% Longitudinal hyperfine; central and common nuclear Zeeman phases refocus
distance=vecnorm(xyz,2,2);
couplings=hyperfine*(1-3*((xyz*field_axis)./distance).^2)./distance.^3;

% Multiply all irreducible pair echoes; each singleton Hahn echo is unity
log_echo=zeros(1,numel(times));
n_spins=size(xyz,1);
pair_tol=abs_tol/max(1,n_spins*(n_spins-1)/2);
max_angle=max(times)/8;
for n=1:(n_spins-1)

    % Nuclear secular Dzz and hyperfine difference for all partners of site n
    delta_r=xyz(n+1:end,:)-xyz(n,:);
    distance=vecnorm(delta_r,2,2);
    pair=dipolar*(1-3*((delta_r*field_axis)./distance).^2)./distance.^3;
    detuning=couplings(n)-couplings(n+1:end);

    % In the flip-flop subspace H_±=(-Dzz σx ± ΔA σz)/4 up to identity
    freq=hypot(pair,detuning);
    active=freq>0;
    freq=freq(active);
    contrast=(2*(pair(active)./freq).*(detuning(active)./freq)).^2;

    % Omit only losses bounded by the per-pair share of the error budget
    keep=contrast.*min(1,(max_angle*freq).^4)>pair_tol;
    if any(keep)

        % Accumulate pair products in log space without discarding weak retained losses
        losses=contrast(keep).*sin(freq(keep).*(times.'/8)).^4;
        log_echo=log_echo+sum(log1p(-losses),1);

    end

end
echo=exp(log_echo).';

end

% Consistency enforcement
function grumble(xyz,field_axis,hyperfine,dipolar,times,abs_tol)
if ~isnumeric(xyz)||~isreal(xyz)||~ismatrix(xyz)||size(xyz,2)~=3||any(~isfinite(xyz(:)))
    error('xyz must be a finite real N-by-3 array.');
end
if any(all(xyz==0,2))||size(unique(xyz,'rows'),1)~=size(xyz,1)
    error('xyz must contain distinct sites away from the central spin.');
end
if ~isnumeric(field_axis)||~isreal(field_axis)||~isequal(size(field_axis),[3 1])||...
   ~(abs(norm(field_axis)-1)<=4*eps)
    error('field_axis must be a real three-element unit column vector.');
end
if ~isnumeric(hyperfine)||~isreal(hyperfine)||~isscalar(hyperfine)||~isfinite(hyperfine)
    error('hyperfine must be a finite real scalar.');
end
if ~isnumeric(dipolar)||~isreal(dipolar)||~isscalar(dipolar)||~isfinite(dipolar)
    error('dipolar must be a finite real scalar.');
end
if ~isnumeric(times)||~isreal(times)||~iscolumn(times)||isempty(times)||...
   any(~isfinite(times))||any(times<0)
    error('times must be a non-empty finite non-negative column vector.');
end
if ~isnumeric(abs_tol)||~isreal(abs_tol)||~isscalar(abs_tol)||~isfinite(abs_tol)||abs_tol<0
    error('abs_tol must be a finite non-negative real scalar.');
end
end


