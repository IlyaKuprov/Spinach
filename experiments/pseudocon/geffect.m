% Effective g-tensor for the user-specified Kramers 
% doublet, computed as described in
% 
%        http://dx.doi.org/10.1063/1.4793736
%
% Syntax:
%
%            g=geffect(spin_system,states)
%
% Parameters:
%
%    states - the numbers of the states to use
%             (numbered sequentially from the 
%              lowest to the highest energy)
%
% Outputs:
%
%    g      - 3x3 g-tensor matrix in Bohr mag-
%             neton units
%
% i.kuprov@soton.ac.uk
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=geffect.m>

function g=geffect(spin_system,states)

% Check consistency
grumble(spin_system,states);

% Get the g-tensor for each spin
for n=spin_system.comp.nspins:-1:1
    g{n}=gtensorof(spin_system,n);   
end

% Get Sx, Sy, Sz operators for each spin
for n=spin_system.comp.nspins:-1:1
    Sx{n}=(operator(spin_system,{'L+'},{n})+...
           operator(spin_system,{'L-'},{n}))/2;
    Sy{n}=(operator(spin_system,{'L+'},{n})-...
           operator(spin_system,{'L-'},{n}))/2i;
    Sz{n}=operator(spin_system,{'Lz'},{n});
end

% Get magnetic moment operators
Mx=sparse(0); My=sparse(0); Mz=sparse(0);
for n=1:spin_system.comp.nspins
    Mx=Mx-g{n}(1,1)*Sx{n}-g{n}(1,2)*Sy{n}-g{n}(1,3)*Sz{n};
    My=My-g{n}(2,1)*Sx{n}-g{n}(2,2)*Sy{n}-g{n}(2,3)*Sz{n};
    Mz=Mz-g{n}(3,1)*Sx{n}-g{n}(3,2)*Sy{n}-g{n}(3,3)*Sz{n};
end

% Get the complete Hamiltonian
[H,Q]=hamiltonian(assume(spin_system,'labframe')); 
H=H+orientation(Q,[0 0 0]); H=(H+H')/2;

% Find the eigenstates
[V,D]=eig(full(H)); D=diag(D);

% Sort the eigenstates
[~,index]=sort(D,'ascend'); V=V(:,index); 

% Pick out the states
V=V(:,states);

% Equation 62 in http://dx.doi.org/10.1063/1.4793736
phi{1}=V'*Mx*V; phi{2}=V'*My*V; phi{3}=V'*Mz*V;

% Equation 61 in http://dx.doi.org/10.1063/1.4793736
for n=1:3
    for k=1:3
        G(n,k)=(phi{n}(1,1)-phi{n}(2,2))*(phi{k}(1,1)-phi{k}(2,2))+...
              2*phi{n}(1,2)*phi{k}(2,1)+2*phi{k}(1,2)*phi{n}(2,1); %#ok<AGROW>
    end
end

% Matrix square root
g=real(sqrtm(G));

end

% Consistency enforcement
function grumble(spin_system,states)
if ~strcmp(spin_system.bas.formalism,'zeeman-hilb')
    error('this function is only available in Hilbert space.');
end
if (~isnumeric(states))||(~isreal(states))||...
   (numel(states)~=2)||any(mod(states,1)~=0)||any(states<1)
    error('the two elements of states argument must be positive integers.');
end
if any(states>size(spin_system.bas.basis,1))
    error('the requested state number exceeds the number of states in the system.');
end
if states(1)==states(2)
    error('the two states must be distinct.');
end
end

% With regard to the idea of whether you have a right to health care, 
% you have to realize what that implies. It's not an abstraction. I'm 
% a physician. That means you have a right to come to my house and 
% conscript me. It means you believe in slavery. It means that you're
% going to enslave not only me, but the janitor at my hospital, the 
% person who cleans my office, the assistants who work in my office,
% the nurses.
%
% Basically, once you imply a belief in a right to someone's services,
% do you have a right to plumbing? Do you have a right to water? Do 
% you have right to food? You're basically saying you believe in sla-
% very. You're saying you believe in taking and extracting from ano-
% ther person. Our founding documents were very clear about this. You
% have a right to pursue happiness but there's no guarantee of physi-
% cal comfort. There's no guarantee of concrete items. In order to 
% give something concrete, you have to take it from someone. So the-
% re's an implied threat of force.
% 
% If I'm a physician in your community and you say you have a right to
% health care, do you have a right to beat down my door with the poli-
% ce, escort me away and force me to take care of you? That's ultima-
% tely what the right to free health care would be. If you believe in 
% a right to health care, you're believing in basically the use of for-
% ce to conscript someone to do your bidding.
%
% Senator Rand Paul

