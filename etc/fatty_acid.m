% Spin system approximating that of a fatty acid. Syntax:
%
%             [sys,inter]=fatty_acid(nprotons)
%
% Parameters:
%
%   nprotons - the number of protons that the spin 
%              system should have
%
% Outputs:
%
%   sys, inter - input data structures for Spinach
%
% ilya.kuprov@weizmann.ac.il
% a.j.allami@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=fatty_acid.m>

function [sys,inter]=fatty_acid(nprotons)

% Check consistency
grumble(nprotons);

% Isotope list
sys.isotopes=cell(1,nprotons);
sys.isotopes(:)={'1H'};

% Chemical shifts
inter.zeeman.scalar=cell(1,nprotons);
inter.zeeman.scalar(1:3)={1.0};                                        % Methyl group +
inter.zeeman.scalar(4:2:end)=num2cell(linspace(1.5,5,(nprotons-3)/2)); % linear ramp along
inter.zeeman.scalar(5:2:end)=num2cell(linspace(1.5,5,(nprotons-3)/2)); % the fatty acid chain

% J-couplings
inter.coupling.scalar=cell(nprotons,nprotons);
inter.coupling.scalar{1,2}=20;
inter.coupling.scalar{2,3}=20;
inter.coupling.scalar{3,1}=20;
inter.coupling.scalar{1,4}=7;
inter.coupling.scalar{1,5}=7;
inter.coupling.scalar{2,4}=7;
inter.coupling.scalar{2,5}=7;
inter.coupling.scalar{3,4}=7;
inter.coupling.scalar{3,5}=7;
for n=4:2:(nprotons-3)
    inter.coupling.scalar{n,n+1}=15+randn();
    inter.coupling.scalar{n,n+2}=7 +randn();
    inter.coupling.scalar{n,n+3}=7 +randn();
end
for n=5:2:(nprotons-2)
    inter.coupling.scalar{n,n+1}=7 +randn();
    inter.coupling.scalar{n,n+2}=7 +randn();
end

end

% Consistency enforcement
function grumble(nprotons)
if (~isnumeric(nprotons))||(~isreal(nprotons))||...
   (mod(nprotons,1)~=0)||(mod((nprotons-3)/2,1)~=0)||...
   ((nprotons-3)/2<1)
    error('(nprotons-3)/2 must be a positive integer.');
end
end

% Sometimes you see beautiful people with no brains. Sometimes
% you have ugly people who are intelligent, like scientists.
%
% Jose Mourinho

