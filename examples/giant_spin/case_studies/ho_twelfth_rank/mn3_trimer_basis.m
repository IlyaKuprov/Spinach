% Effective basis of the (CH6N3)2MnCl4 trimer case studies, built as in
% Section 3.2 of https://arxiv.org/abs/2609.16352: the common eigenstates
% of the isotropic exchange Hamiltonian and the total S_z, obtained by
% diagonalising the exchange Hamiltonian with a perturbative 1 mT Zeeman
% term, are selected either as the lowest state of every total S_z pro-
% jection (16 states, Table 1 of the paper) or as every state within
% 25 cm^-1 of the ground state (26 states, Table 2 of the paper). Syntax:
%
%                [P,msz]=mn3_trimer_basis(spin_system,nstates)
%
% Parameters:
%
%    spin_system - Spinach spin system of the trimer, three E6 spins
%                  in zeeman-hilb formalism
%
%    nstates     - 16 or 26, the basis set of the paper to return
%
% Outputs:
%
%    P           - matrix with the basis states in columns, 216 x nstates
%
%    msz         - column of total S_z projections of the basis states
%
% ilya.kuprov@weizmann.ac.il

function [P,msz]=mn3_trimer_basis(spin_system,nstates)

% Check consistency
grumble(spin_system,nstates);

% Isotropic exchange Hamiltonian and
% the Zeeman operator per Tesla
I=hamiltonian(assume(spin_system,'labframe'));
Z=hamiltonian(assume(spin_system,'labframe','zeeman'));
H_exch=I-Z;

% Common eigenstates of the exchange Hamiltonian
% and the total S_z from a 1 mT Zeeman term
H=full(H_exch+1e-3*Z); [V,E]=eig((H+H')/2,'vector');
Sz=full(operator(spin_system,'Lz','E6'));
msz=round(2*real(diag(V'*Sz*V)))/2;

% Energies relative to the ground state, cm^-1
E=hz2icm((E-min(E))/(2*pi));

% Lowest state of every projection, or every
% state within 25 cm^-1 of the ground state
if nstates==16
    sel=zeros(16,1); proj=(-7.5:7.5)';
    for n=1:16
        idx=find(msz==proj(n)); [~,k]=min(E(idx)); sel(n)=idx(k);
    end
else
    sel=find(E<25);
end
P=V(:,sel); msz=msz(sel);

end

% Consistency enforcement
function grumble(spin_system,nstates)
if ~strcmp(spin_system.bas.formalism,'zeeman-hilb')
    error('this function is only available in zeeman-hilb formalism.');
end
if ~isequal(spin_system.comp.isotopes,{'E6','E6','E6'})
    error('spin_system must contain three E6 spins.');
end
if (~isnumeric(nstates))||(~isscalar(nstates))||(~ismember(nstates,[16 26]))
    error('nstates must be 16 or 26.');
end
end

