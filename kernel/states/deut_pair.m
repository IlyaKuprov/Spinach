% All possible states of a spin-1 pair, classified by the total spin 
% into singlet, triplet, and quartet. Syntax:
%
%        [S,T,Q,ODS,PDS]=deut_pair(spin_system,spin_a,spin_b)
%
% Arguments:
%
%   spin_a  - the number of the first spin 
%
%   spin_b  - the number of the second spin 
%
%   options.dephasing - set to 1 to eliminate the states that are
%                       not stationary under Az+Bz Hamiltonian, 
%                       the default is to keep everything
%
% Outputs:
%
%   S       - singet state density matrix (Hilbert space)
%             or state vector (Liouville space)
%
%   T       - triplet state density matrices (Hilbert space)
%             or state vectors (Liouville space), ordered in
%             a cell array as {T+,T0,T-}
%
%   Q       - quintet state density matrices (Hilbert space)
%             or state vectors (Liouville space), ordered in
%             a cell array as {Q++,Q+,Q0,Q-,Q--}
%
% WARNING: the states above are NOT irreducible spherical tensors - 
%          Bargon just Kroneckered up some Zeeman states and gave
%          them what looked to him like appropriate labels.
%
% i.kuprov@soton.ac.uk
% anakin.aden@mpinat.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=deut_pair.m>

function [S,T,Q]=deut_pair(spin_system,spin_a,spin_b,options)

% Set defaults
if ~exist('options','var'), options.dephasing=0; end

% Check consistency
grumble(spin_system,spin_a,spin_b,options);

% Build component vectors in Hilbert space as per Eq 1
% in https://doi.org/10.1016/S0009-2614(98)00784-2
alp=[1; 0; 0]; bet=[0; 1; 0]; gam=[0; 0; 1];
S0=(1/sqrt(3))*(kron(alp,gam)-kron(bet,bet)+kron(gam,alp));
Tp=(1/sqrt(2))*(kron(alp,bet)-kron(bet,alp));
T0=(1/sqrt(2))*(kron(alp,gam)-kron(gam,alp));
Tm=(1/sqrt(2))*(kron(bet,gam)-kron(gam,bet));
Qpp=kron(alp,alp); 
Qp=(1/sqrt(2))*(kron(alp,bet)+kron(bet,alp));
Q0=(1/sqrt(6))*(kron(alp,gam)+2*kron(bet,bet)+kron(gam,alp));
Qm=(1/sqrt(2))*(kron(bet,gam)+kron(gam,bet));
Qmm=kron(gam,gam);

% Obtain spherical tensor expansions
S=0; T={0,0,0}; Q={0,0,0,0,0}; IST=irr_sph_ten(3);
for n=1:numel(IST)
    for k=1:numel(IST)

        % Get spherical tensor indices
        [L1,M1]=lin2lm(n-1); [L2,M2]=lin2lm(k-1);

        % Skip non-stationary states on user request
        if (M1~=0)&&(M2~=0)&&(options.dephasing==1), continue; end

        % Build Spinach state descriptor
        descr_a=['T' int2str(L1) ',' int2str(M1)];
        descr_b=['T' int2str(L2) ',' int2str(M2)];

        % Get the Spinach state
        rho=state(spin_system,{descr_a,descr_b},...
                              {spin_a, spin_b });

        % Get and normalise the two-spin state
        tss=kron(IST{n},IST{k}); tss=tss/norm(tss,'fro');

        % Normalise the Spinach state
        switch spin_system.bas.formalism

            case 'zeeman-hilb'

                % Hilbert space
                rho=rho/norm(rho,'fro');
                
            case {'zeeman-liouv','sphten-liouv'}

                % Liouville space
                rho=rho/norm(full(rho),2);

            otherwise

                % Complain and bomb out
                error('unsupported formalism.');

        end

        % Add to the relevant totals
        S   =S   +trace(S0' *tss'*S0 )*rho;
        T{1}=T{1}+trace(Tp' *tss'*Tp )*rho;
        T{2}=T{2}+trace(T0' *tss'*T0 )*rho;
        T{3}=T{3}+trace(Tm' *tss'*Tm )*rho;
        Q{1}=Q{1}+trace(Qpp'*tss'*Qpp)*rho;
        Q{2}=Q{2}+trace(Qp' *tss'*Qp )*rho;
        Q{3}=Q{3}+trace(Q0' *tss'*Q0 )*rho;
        Q{4}=Q{4}+trace(Qm' *tss'*Qm )*rho;
        Q{5}=Q{5}+trace(Qmm'*tss'*Qmm)*rho;

    end

end

end

% Consistency enforcement
function grumble(spin_system,spin_a,spin_b,options)
if (~isnumeric(spin_a))||(~isnumeric(spin_b))||...
   (~isscalar(spin_a))||(~isscalar(spin_b))||...
   (mod(spin_a,1)~=0)||(mod(spin_b,1)~=0)||...
   (spin_a<1)||(spin_b<1)||(spin_a==spin_b)
    error('spin indices must be different positive integers.');
end
if (spin_a>spin_system.comp.nspins)||...
   (spin_b>spin_system.comp.nspins)
    error('spin index exceeds the number of spins.');
end
if (spin_system.comp.mults(spin_a)~=3)||...
   (spin_system.comp.mults(spin_b)~=3)
    error('both particles must have spin 1.')
end
if ~isfield(options,'dephasing')
    error('options.dephasing must be set.');
end
if (~isnumeric(options.dephasing))||...
   (~isscalar(options.dephasing))||...
   (~ismember(options.dephasing,[0 1]))
    error('options.dephasing must be 1 or 0.');
end
end

% My suggestion was quite simple: put that needed code number in a
% little capsule, and then implant that capsule right next to the
% heart of a volunteer. The volunteer would carry with him a big,
% heavy butcher knife as he accompanies the President. If ever the
% President wanted to fire nuclear weapons, the only way he could
% do so would be for him first, with his own hands, to kill one hu-
% man being. The President says, "George, I'm sorry but tens of mil-
% lions must die." He has to look at someone and realize what death
% is - what an innocent death is. Blood on the White House carpet. 
% Its reality brought home.
%
% When I suggested this to friends in the Pentagon they said, "My
% God, that's terrible. Having to kill someone would distort the
% President's judgment. He might never push the button."
%
% Roger Fisher

