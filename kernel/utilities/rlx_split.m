% Splits a relaxation superoperator into longitudinal, trans-
% verse and mixed components. Syntax:
%
%             [R1,R2,Rm]=rlx_split(spin_system,R)
%
% Parameters:
%
%      R   - a relaxation superoperator in sphten-liouv
%            formalism
%
% Outputs:
%
%      R1  - the part of R acting on purely longitudinal
%            single-spin states
%
%      R2  - the part of R acting on purely transverse
%            single-spin states
%
%      Rm  - the rest of R
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rlx_split.m>

function [R1,R2,Rm]=rlx_split(spin_system,R)

% Check consistency
grumble(spin_system,R);

% Interpret the basis
[L,M]=lin2lm(spin_system.bas.basis);

% Index single- and multi-spin orders (sso and mso)
sso_mask=(sum(logical(spin_system.bas.basis),2)==1);
mso_mask=(sum(logical(spin_system.bas.basis),2)>1 );

% Index longlitudinal and transverse states
long_sso_mask=any((L>0)&(M==0),2)&sso_mask;
tran_sso_mask=any((L>0)&(M~=0),2)&sso_mask;

% Split the relaxation superoperator
R1=R; R1(~long_sso_mask,~long_sso_mask)=0;
R2=R; R2(~tran_sso_mask,~tran_sso_mask)=0;
Rm=R; R2(~mso_mask,~mso_mask)=0;

end

% Consistency enforcement
function grumble(spin_system,R)
if (~isfield(spin_system,'bas'))||(~isfield(spin_system.bas,'formalism'))
    error('spin_system object does not contain the necessary information.');
end
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('this function is only available in sphten-liouv formalism.');
end
if (~isnumeric(R))||(~ismatrix(R))||(size(R,1)~=size(R,2))
    error('R must be a square matrix.');
end
end

% He that reproveth a scorner getteth to himself shame: and he that
% rebuketh a wicked man getteth himself a blot. Reprove not a scorn-
% er, lest he hate thee: rebuke a wise man, and he will love thee.
% Give instruction to a wise man, and he will be yet wiser: teach a
% just man, and he will increase in learning.
%
% Proverbs 9:7

