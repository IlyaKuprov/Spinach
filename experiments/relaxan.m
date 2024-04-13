% Automated relaxation theory analysis. Prints longitudinal
% and transverse relaxation rates and times for all spins in
% the system. Syntax:
%
%     [r1,r2,t1,t2,R]=relaxan(spin_system,euler_angles)
%
% Parameters:
%
%      euler_angles - optional euler angles for situations
%                     when relaxation properties are orien-
%                     tation-dependent
%
% Outputs:
%
%      r1 - a vector of longitudinal relaxation rates 
%           for each spin
%
%      r2 - a vector of transverse relaxation rates 
%           for each spin
%
%      t1 - a vector of longitudinal relaxation times 
%           for each spin
%
%      t2 - a vector of transverse relaxation times 
%           for each spin
%
%       R - complete relaxation superoperator
%
% Note: dynamic frequency shifts are dropped.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=relaxan.m>

function [r1,r2,t1,t2,R]=relaxan(spin_system,euler_angles)

% Check consistency
grumble(spin_system);

% Compute the relaxation superoperator
if exist('euler_angles','var')
    R=relaxation(spin_system,euler_angles);
else
    R=relaxation(spin_system);
end

% Preallocate outputs
r1=zeros(spin_system.comp.nspins,1);
r2=zeros(spin_system.comp.nspins,1);
t1=zeros(spin_system.comp.nspins,1);
t2=zeros(spin_system.comp.nspins,1);

% Fill the arrays
parfor n=1:spin_system.comp.nspins
    Lz=state(spin_system,{'Lz'},{n},'cheap');
    r1(n)=-real((Lz'*R*Lz)/(Lz'*Lz)); t1(n)=1/r1(n);
    Lp=state(spin_system,{'L+'},{n},'cheap');
    r2(n)=-real((Lp'*R*Lp)/(Lp'*Lp)); t2(n)=1/r2(n);
end

% Do the printing
report(spin_system,'===============================================================');
report(spin_system,' Number  Isotope  R1(Hz)      R2(Hz)      T1(s)       T2(s)    ');
report(spin_system,'===============================================================');
for n=1:spin_system.comp.nspins
    report(spin_system,[' ' pad(num2str(n),8) pad(spin_system.comp.isotopes{n},9)...
                        pad(num2str(r1(n),'%4.3e'),12) pad(num2str(r2(n),'%4.3e'),12)...
                        pad(num2str(t1(n),'%4.3e'),12) pad(num2str(t2(n),'%4.3e'),12)...
                        spin_system.comp.labels{n}]);
end
report(spin_system,'===============================================================');

end

% Consistency enforcement
function grumble(spin_system)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function is only available in Liouville space.');
end
end

% By the grace of reality and the nature of life, man - every man - is an
% end in himself, he exists for his own sake, and the achievement of his
% own happiness is his highest moral purpose.
%
% Ayn Rand, "Atlas Shrugged"

