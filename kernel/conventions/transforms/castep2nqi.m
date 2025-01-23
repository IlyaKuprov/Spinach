% Converts CASTEP EFG tensor (it is printed in atomic units) to NQI
% 3x3 tensor in Hz that is required by Spinach. Syntax:
%
%                     nqi=castep2nqi(V,Q,I)
%
% Parameters:
%
%     V   - EFG tensor from CASTEP output, a.u.
%
%     Q   - nuclear quadrupole moment, barn
%
%     I   - nuclear spin quantum number
%
% Outputs:
%
%    nqi  - 3x3 matrix in Hz, ready for input 
%           into create.m function
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=castep2nqi.m>

function nqi=castep2nqi(V,Q,I)

% Check consistency
grumble(V,Q,I)

% Fundamental constants
efg_atomic=9.717362e+21;
e_charge=1.60217657e-19;
h_planck=6.62606957e-34;

% Calculation
nqi=V*efg_atomic*(Q*1e-28)*e_charge/(h_planck*2*I*(2*I-1));

end

% Consistency enforcement
function grumble(V,Q,I)
if (~isnumeric(V))||(~isnumeric(Q))||(~isnumeric(I))
    error('all inputs must be numeric.');
end
if (~isreal(V))||(~isreal(Q))||(~isreal(I))
    error('all inputs must be real.');
end
if ~all(size(V)==[3 3])
    error('V argument must be a 3x3 matrix.');
end
if numel(Q)~=1
    error('Q argument must have a single element.');
end
if (numel(I)~=1)||(I<1)||(mod(2*I+1,1)~=0)
    error('I must be an integer or half-integer greater or equal to 1.');
end
end

% My children! We have fought in many battles together, over mountaintops
% and beach heads, through forests and deserts. I have seen great acts of
% valor from each one of you, which does my heart proud. I have also seen
% dirty fighting, backstabbing, cruel and wanton feats of savagery, which
% pleases me equally well. For you are all warriors.
%
% Queen Potema

