% Isotropic second-order quadrupolar shift of a Zeeman transition.
% Syntax:
%
%  shift=quad_iso_shift(I,cq_hz,eta,nu0_hz,m_from,m_to)
%
% Parameters:
%
%  I             - integer or half-integer nuclear spin, at least one
%  cq_hz         - quadrupolar coupling e^2*q*Q/h in Hz
%  eta           - asymmetry, between zero and one
%  nu0_hz        - signed coefficient of Iz in the Zeeman Hamiltonian in Hz
%  m_from,m_to   - distinct allowed magnetic quantum numbers
%
% Outputs:
%
%  shift         - second-order correction to E(m_to)-E(m_from), in Hz
%
% The orientationally averaged squared off-diagonal quadrupolar matrix
% elements share the factor trace(Q^2)/20, where Q is the interaction tensor
% in Hz. Perturbation denominators use nu0_hz*(m-m_prime).
%
% talos@spindynamics.org

function shift=quad_iso_shift(I,cq_hz,eta,nu0_hz,m_from,m_to)

% Check consistency
grumble(I,cq_hz,eta,nu0_hz,m_from,m_to);

% Form upward and downward squared ladder factors from each level
m=[m_from m_to];
up=(I-m).*(I+m+1);
down=(I+m).*(I-m+1);
up_next=(I-m-1).*(I+m+2);
down_next=(I+m-1).*(I-m+2);

% Sum the denominators for the two single- and two double-quantum paths
energy=down.*(2*m-1).^2-up.*(2*m+1).^2+...
       (down.*down_next-up.*up_next)/2;

% Apply the uniform orientation average and subtract the two levels
tensor_zz=cq_hz/(2*I*(2*I-1));
shift=tensor_zz^2*(3+eta^2)*diff(energy)/(40*nu0_hz);

end

% Consistency enforcement
function grumble(I,cq_hz,eta,nu0_hz,m_from,m_to)
if ~isnumeric(I)||~isreal(I)||~isscalar(I)||I<1||mod(2*I,1)~=0
    error('I must be an integer or half-integer at least one.');
end
if ~isnumeric(cq_hz)||~isreal(cq_hz)||~isscalar(cq_hz)||~isfinite(cq_hz)
    error('cq_hz must be finite and real.');
end
if ~isnumeric(eta)||~isreal(eta)||~isscalar(eta)||~isfinite(eta)||eta<0||eta>1
    error('eta must be between zero and one.');
end
if ~isnumeric(nu0_hz)||~isreal(nu0_hz)||~isscalar(nu0_hz)||~isfinite(nu0_hz)||nu0_hz==0
    error('nu0_hz must be finite, nonzero, and real.');
end
if ~isnumeric(m_from)||~isnumeric(m_to)||~isreal(m_from)||~isreal(m_to)||...
   ~isscalar(m_from)||~isscalar(m_to)||m_from==m_to||...
   any(abs([m_from m_to])>I)||any(mod(I-[m_from m_to],1)~=0)
    error('m_from and m_to must be distinct Zeeman levels of I.');
end
end


