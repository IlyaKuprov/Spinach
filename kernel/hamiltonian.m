% Hamiltonian operator or superoperator and its rotational decomposi-
% tion. Descriptor and operator generation are parallelised. Syntax:
%
%             [H,Q]=hamiltonian(spin_system,operator_type)
%
% Parameters: 
%
%    in Liouville space, operator_type can be set to
%                      
%          'left' - produces left side product superoperator
%
%         'right' - produces right side product superoperator
%
%          'comm' - produces commutation superoperator (default)
%
%         'acomm' - produces anticommutation superoperator
%
%    in Hilbert space this parameter is ignored.
%
% Outputs:
%
%     H   - rotationally invariant part of the Hamiltonian
%
%     Q   - irreducible components of the anisotropic part,
%           use orientation.m to get the full Hamiltonian
%           at each specific orientation
%
% Note: the code has a few rather eccentric blocks that bring the 
%       memory footprint to the absolute minimum and work around 
%       the sparse matrix addition efficiency problem.
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
% hannah.hogben@chem.ox.ac.uk
% d.savostyanov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=hamiltonian.m>

function [H,Q]=hamiltonian(spin_system,operator_type)

% Set the default for the type
if ~exist('operator_type','var'), operator_type='comm'; end

% Check consistency
grumble(spin_system,operator_type);

% Decide if the Q part is required
build_aniso=(nargout>1);

% Inform the user
report(spin_system,'building Hamiltonian descriptor...');

% Preallocate spin number tables
nL=zeros(spin_system.comp.nspins,8);
nS=zeros(spin_system.comp.nspins,8);

% Preallocate operator specifications
opL(1:spin_system.comp.nspins,1:8)={'E'};
opS(1:spin_system.comp.nspins,1:8)={'E'};

% Preallocate isotropic Hamiltonian coefficients
isotropic=zeros(spin_system.comp.nspins,8);

% Preallocate spherical tensor coefficients
ist_coeff{1}=zeros([spin_system.comp.nspins 8 3],'like',1i); % Rank 1 [+1, 0, -1]
ist_coeff{2}=zeros([spin_system.comp.nspins 8 5],'like',1i); % Rank 2 [+2, +1, 0, -1, 2]

% Preallocate irreducible components
irr_comp{1}=zeros([spin_system.comp.nspins 8 3],'like',1i);  % Rank 1 [+1, 0, -1]
irr_comp{2}=zeros([spin_system.comp.nspins 8 5],'like',1i);  % Rank 2 [+2, +1, 0, -1, 2]

% Process Zeeman interactions and NQI
for n=1:spin_system.comp.nspins
    
    % Write the isotropic Zeeman part
    switch spin_system.inter.zeeman.strength{n}
        
        case 'full'
            
            % Keep the carrier frequency
            zeeman_iso=trace(spin_system.inter.zeeman.matrix{n})/3;
            
            % Update the Hamiltonian
            if abs(zeeman_iso)>spin_system.tols.liouv_zero
                
                % Inform the user
                report(spin_system,['complete isotropic Zeeman interaction for spin ' num2str(n) '...']);
                report(spin_system,['           (Lz) x ' num2str(zeeman_iso/(2*pi)) ' Hz']);
                
                % Update the Hamiltonian descriptor
                nL(n,2)=n; opL(n,2)={'Lz'}; isotropic(n,2)=zeeman_iso;
                
            end
            
        case 'secular'
            
            % Subtract the carrier frequency
            zeeman_iso=trace(spin_system.inter.zeeman.matrix{n})/3-...
                             spin_system.inter.basefrqs(n);
            
            % Update the Hamiltonian
            if abs(zeeman_iso)>spin_system.tols.liouv_zero
                
                % Inform the user
                report(spin_system,['offset isotropic Zeeman interaction for spin ' num2str(n) '...']);
                report(spin_system,['           (Lz) x ' num2str(zeeman_iso/(2*pi)) ' Hz']);
                
                % Update the Hamiltonian descriptor
                nL(n,2)=n; opL(n,2)={'Lz'}; isotropic(n,2)=zeeman_iso;
                
            end
            
        case 'ignore'
            
            % Inform the user
            report(spin_system,['isotropic Zeeman interaction ignored for spin ' num2str(n) '.']);
            
        otherwise
            
            % Bomb out with unexpected strength parameters
            error(['unknown strength specification for the Zeeman interaction of spin ' num2str(nL)]);
            
    end
    
    % Process anisotropic Zeeman and NQI if required
    if build_aniso
        
        % Get spherical tensor components for Zeeman interactions
        [~,lam_zeeman,phi_zeeman]=mat2sphten(spin_system.inter.zeeman.matrix{n});
        
        % Only process significant Zeeman anisotropies
        if (norm(lam_zeeman,2)>spin_system.tols.liouv_zero)||...
           (norm(phi_zeeman,2)>spin_system.tols.liouv_zero)
            
            % Store coefficients
            for k=1:3
                irr_comp{1}(n,k,:)=lam_zeeman;
                irr_comp{2}(n,k,:)=phi_zeeman;
            end
            
            % Write irreducible spherical tensors
            switch spin_system.inter.zeeman.strength{n}
                
                case 'full'
                    
                    % Inform the user
                    report(spin_system,['complete anisotropic Zeeman interaction for spin ' num2str(n) '...']);
                    if norm(lam_zeeman,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[1,+1]:      -0.5*(Lp) x ' num2str(lam_zeeman(1)/(2*pi))   ' Hz']);
                        report(spin_system,['  T[1,-1]:      -0.5*(Lm) x ' num2str(lam_zeeman(3)/(2*pi))   ' Hz']);
                    end
                    if norm(phi_zeeman,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[2,+1]:      -0.5*(Lp) x ' num2str(phi_zeeman(2)/(2*pi))   ' Hz']);
                        report(spin_system,['  T[2, 0]: sqrt(2/3)*(Lz) x ' num2str(phi_zeeman(3)/(2*pi))   ' Hz']);
                        report(spin_system,['  T[2,-1]:       0.5*(Lm) x ' num2str(phi_zeeman(4)/(2*pi))   ' Hz']);
                    end
                    
                    % Prepare spherical tensor descriptors
                    nL(n,1)=n; opL(n,1)={'L+'}; ist_coeff{1}(n,1,1)=-0.5; ist_coeff{2}(n,1,2)=-0.5;
                    nL(n,2)=n; opL(n,2)={'Lz'};                           ist_coeff{2}(n,2,3)=sqrt(2/3);
                    nL(n,3)=n; opL(n,3)={'L-'}; ist_coeff{1}(n,3,3)=-0.5; ist_coeff{2}(n,3,4)=+0.5;
                    
                case 'secular'
                    
                    % Inform the user
                    report(spin_system,['Z part of the anisotropic Zeeman interaction for spin ' num2str(n) '...']);
                    report(spin_system,['  T[2, 0]: sqrt(2/3)*(Lz) x ' num2str(phi_zeeman(3)/(2*pi)) ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    nL(n,2)=n; opL(n,2)={'Lz'};                           ist_coeff{2}(n,2,3)=sqrt(2/3);
                    
                case 'ignore'
                    
                    % Inform the user
                    report(spin_system,['anisotropic Zeeman interaction ignored for spin ' num2str(n) '.']);
                    
                otherwise
                    
                    % Bomb out with unexpected strength parameters
                    error(['unknown Zeeman interaction strength specification for spin ' num2str(n) '.']);
                    
            end
            
        end
        
        % Process quadrupolar interactions
        if norm(spin_system.inter.coupling.matrix{n,n},2)>2*pi*spin_system.tols.inter_cutoff
            
            % Get spherical tensor components for NQIs
            [iso_quad,lam_quad,phi_quad]=mat2sphten(spin_system.inter.coupling.matrix{n,n});
            
            % Catch disallowed components
            if (norm(iso_quad,2)>1e-6*norm(phi_quad,2))||...
               (norm(lam_quad,2)>1e-6*norm(phi_quad,2))
                error('quadratic couplings must be traceless and symmetric.');
            end
            
            % Store coefficients
            for k=4:8
                irr_comp{1}(n,k,:)=[0 0 0];
                irr_comp{2}(n,k,:)=phi_quad;
            end
            
            % Process the coupling
            switch spin_system.inter.coupling.strength{n,n}
                
                case 'strong'
                    
                    % Inform the user
                    report(spin_system,['complete quadratic coupling for spin ' num2str(n) '...']);
                    report(spin_system,['           T(2,+2) x ' num2str(phi_quad(1)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['           T(2,+1) x ' num2str(phi_quad(2)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['           T(2, 0) x ' num2str(phi_quad(3)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['           T(2,-1) x ' num2str(phi_quad(4)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['           T(2,-2) x ' num2str(phi_quad(5)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    nL(n,4)=n; opL(n,4)={'T2,+2'}; ist_coeff{2}(n,4,1)=1;
                    nL(n,5)=n; opL(n,5)={'T2,+1'}; ist_coeff{2}(n,5,2)=1;
                    nL(n,6)=n; opL(n,6)={'T2,0'};  ist_coeff{2}(n,6,3)=1;
                    nL(n,7)=n; opL(n,7)={'T2,-1'}; ist_coeff{2}(n,7,4)=1;
                    nL(n,8)=n; opL(n,8)={'T2,-2'}; ist_coeff{2}(n,8,5)=1;
                    
                case 'secular'
                    
                    % Inform the user
                    report(spin_system,['secular part of the quadratic coupling for spin ' num2str(n) '...']);
                    report(spin_system,['           T(2, 0) x ' num2str(phi_quad(3)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    nL(n,6)=n; opL(n,6)={'T2,0'}; ist_coeff{2}(n,6,3)=1;
                    
                case 'ignore'
                    
                    % Inform the user
                    report(spin_system,['quadratic coupling ignored for spin ' num2str(n) '.']);
                    
                otherwise
                    
                    % Bomb out with unexpected strength parameters
                    error(['unknown strength specification for the quadratic coupling of spin ' num2str(n)]);
                    
            end
            
        end
        
    end
    
end

% Pack single-spin descriptor table
D1=table(reshape(nL, [8*spin_system.comp.nspins 1]),...
         reshape(nS, [8*spin_system.comp.nspins 1]),...
         reshape(opL,[8*spin_system.comp.nspins 1]),...
         reshape(opS,[8*spin_system.comp.nspins 1]),...
         reshape(isotropic,[8*spin_system.comp.nspins 1]),...
         [reshape(ist_coeff{1},[8*spin_system.comp.nspins 3])...
          reshape(ist_coeff{2},[8*spin_system.comp.nspins 5])],...
         [reshape(irr_comp{1}, [8*spin_system.comp.nspins 3])...
          reshape(irr_comp{2}, [8*spin_system.comp.nspins 5])],...
         'VariableNames',{'nL','nS','opL','opS','isotropic',...
                          'ist_coeff','irr_comp'});

% Kill insignificant rows
mask_a=(abs(D1.isotropic)       <spin_system.tols.liouv_zero);
mask_b=(sum(abs(D1.irr_comp),2) <spin_system.tols.liouv_zero);
mask_c=(sum(abs(D1.ist_coeff),2)<spin_system.tols.liouv_zero);
D1(mask_a&(mask_b|mask_c),:)=[];

% Tightest possible variable clean-up
clear('nL','nS','opL','opS','isotropic','ist_coeff','irr_comp',...
      'mask_a','mask_b','mask_c','zeeman_iso','n');

% Discover significant coupling tensors and build interacting pair list
[L,S]=find(cellfun(@(x)norm(x,2),spin_system.inter.coupling.matrix)>...
                                 2*pi*spin_system.tols.inter_cutoff);
quad_couplings=(L==S); L(quad_couplings)=[]; S(quad_couplings)=[];
pair_list=[L S]; if isempty(pair_list), pair_list=[]; end

% Preallocate spin indices
nL=zeros(size(pair_list,1),3,3);
nS=zeros(size(pair_list,1),3,3);

% Preallocate operator specifications
opL(1:size(pair_list,1),1:3,1:3)={'E'};
opS(1:size(pair_list,1),1:3,1:3)={'E'};

% Preallocate isotropic Hamiltonian coefficients
isotropic=zeros([size(pair_list,1) 3 3],'like',1i);

% Preallocate spherical tensor coefficients
ist_coeff{1}=zeros(size(pair_list,1),3,3,3);
ist_coeff{2}=zeros(size(pair_list,1),3,3,5);

% Preallocate ireducible components
irr_comp{1}=zeros(size(pair_list,1),3,3,3);
irr_comp{2}=zeros(size(pair_list,1),3,3,5);

% Loop over the pair list
for n=1:size(pair_list,1)
    
    % Extract spin numbers
    L=pair_list(n,1); S=pair_list(n,2);
    
    % Get the isotropic coupling constant
    coupling_iso=trace(spin_system.inter.coupling.matrix{L,S})/3;
    
    % Filter out small isotropic parts
    if abs(coupling_iso)>2*pi*spin_system.tols.inter_cutoff
        
        % Process the coupling
        switch spin_system.inter.coupling.strength{L,S}
            
            case {'strong','secular'}
                
                % Inform the user
                report(spin_system,['complete isotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                report(spin_system,['           (LxSx+LySy+LzSz) x ' num2str(coupling_iso/(2*pi)) ' Hz']);
                
                % Update the Hamiltonian descriptor
                nL(n,2,2)=L; nS(n,2,2)=S; opL(n,2,2)={'Lz'}; opS(n,2,2)={'Lz'}; isotropic(n,2,2)=coupling_iso;
                nL(n,1,3)=L; nS(n,1,3)=S; opL(n,1,3)={'L+'}; opS(n,1,3)={'L-'}; isotropic(n,1,3)=coupling_iso/2;
                nL(n,3,1)=L; nS(n,3,1)=S; opL(n,3,1)={'L-'}; opS(n,3,1)={'L+'}; isotropic(n,3,1)=coupling_iso/2;
                
            case {'z*','*z','zz'}
                
                % Inform the user
                report(spin_system,['(z,z) part of the isotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                report(spin_system,['           (LzSz) x ' num2str(coupling_iso/(2*pi)) ' Hz']);
                
                % Update the Hamiltonian descriptor
                nL(n,2,2)=L; nS(n,2,2)=S; opL(n,2,2)={'Lz'}; opS(n,2,2)={'Lz'}; isotropic(n,2,2)=coupling_iso;
                
            case {'ignore','z+','z-','+z','-z','T(L,+1)','T(L,-1)'}
                
                % Inform the user
                report(spin_system,['isotropic coupling ignored for spins ' num2str(L) ',' num2str(S) '.']);
                
            otherwise
                
                % Bomb out with unexpected strength parameters
                error(['unknown strength specification for the bilinear coupling between spins ' num2str(L) ' and ' num2str(S)]);
                
        end
        
    end
    
    % Process anisotropic part if required
    if build_aniso
        
        % Get spherical tensor components
        [~,lam_coupling,phi_coupling]=mat2sphten(spin_system.inter.coupling.matrix{L,S});
        
        % Only process significant interaction anisotropies
        if (norm(lam_coupling,2)>spin_system.tols.liouv_zero)||...
           (norm(phi_coupling,2)>spin_system.tols.liouv_zero)
        
            % Store coefficients
            for k=1:3
                for m=1:3
                    irr_comp{1}(n,k,m,:)=lam_coupling;
                    irr_comp{2}(n,k,m,:)=phi_coupling;
                end
            end
            
            % Write irreducible spherical tensors
            switch spin_system.inter.coupling.strength{L,S}
                
                case 'strong'
                    
                    % Inform the user
                    report(spin_system,['complete anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    if norm(lam_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[1,+1]:                  -0.5*(LpSz-LzSp) x ' num2str(lam_coupling(1)/(2*pi),'%+0.5e') ' Hz']);
                        report(spin_system,['  T[1, 0]:            -sqrt(1/8)*(LpSm-LmSp) x ' num2str(lam_coupling(2)/(2*pi),'%+0.5e') ' Hz']);
                        report(spin_system,['  T[1,-1]:                  -0.5*(LmSz-LzSm) x ' num2str(lam_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    if norm(phi_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[2,+2]:                        0.5*(LpSp) x ' num2str(phi_coupling(1)/(2*pi),'%+0.5e') ' Hz']);
                        report(spin_system,['  T[2,+1]:                  -0.5*(LzSp+LpSz) x ' num2str(phi_coupling(2)/(2*pi),'%+0.5e') ' Hz']);
                        report(spin_system,['  T[2, 0]: sqrt(2/3)*(LzSz-0.25*(LpSm+LmSp)) x ' num2str(phi_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                        report(spin_system,['  T[2,-1]:                   0.5*(LzSm+LmSz) x ' num2str(phi_coupling(4)/(2*pi),'%+0.5e') ' Hz']);
                        report(spin_system,['  T[2,-2]:                        0.5*(LmSm) x ' num2str(phi_coupling(5)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    
                    % Prepare spherical tensor descriptors
                    nL(n,2,2)=L; nS(n,2,2)=S; opL(n,2,2)={'Lz'}; opS(n,2,2)={'Lz'};                                   ist_coeff{2}(n,2,2,3)=+sqrt(2/3);
                    nL(n,1,3)=L; nS(n,1,3)=S; opL(n,1,3)={'L+'}; opS(n,1,3)={'L-'}; ist_coeff{1}(n,1,3,2)=-sqrt(1/8); ist_coeff{2}(n,1,3,3)=-sqrt(2/3)/4;
                    nL(n,3,1)=L; nS(n,3,1)=S; opL(n,3,1)={'L-'}; opS(n,3,1)={'L+'}; ist_coeff{1}(n,3,1,2)=+sqrt(1/8); ist_coeff{2}(n,3,1,3)=-sqrt(2/3)/4;
                    nL(n,2,1)=L; nS(n,2,1)=S; opL(n,2,1)={'Lz'}; opS(n,2,1)={'L+'}; ist_coeff{1}(n,2,1,1)=+1/2;       ist_coeff{2}(n,2,1,2)=-1/2;
                    nL(n,1,2)=L; nS(n,1,2)=S; opL(n,1,2)={'L+'}; opS(n,1,2)={'Lz'}; ist_coeff{1}(n,1,2,1)=-1/2;       ist_coeff{2}(n,1,2,2)=-1/2;
                    nL(n,2,3)=L; nS(n,2,3)=S; opL(n,2,3)={'Lz'}; opS(n,2,3)={'L-'}; ist_coeff{1}(n,2,3,3)=+1/2;       ist_coeff{2}(n,2,3,4)=+1/2;
                    nL(n,3,2)=L; nS(n,3,2)=S; opL(n,3,2)={'L-'}; opS(n,3,2)={'Lz'}; ist_coeff{1}(n,3,2,3)=-1/2;       ist_coeff{2}(n,3,2,4)=+1/2;
                    nL(n,1,1)=L; nS(n,1,1)=S; opL(n,1,1)={'L+'}; opS(n,1,1)={'L+'};                                   ist_coeff{2}(n,1,1,1)=+1/2;
                    nL(n,3,3)=L; nS(n,3,3)=S; opL(n,3,3)={'L-'}; opS(n,3,3)={'L-'};                                   ist_coeff{2}(n,3,3,5)=+1/2;
                    
                case 'z*'
                    
                    % Inform the user
                    report(spin_system,['(z,*) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    if norm(lam_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[1,+1]:                0.5*(LzSp) x ' num2str(lam_coupling(1)/(2*pi),'%+0.5e') ' Hz']);
                        report(spin_system,['  T[1,-1]:                0.5*(LzSm) x ' num2str(lam_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    if norm(phi_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[2,+1]:               -0.5*(LzSp) x ' num2str(phi_coupling(2)/(2*pi),'%+0.5e') ' Hz']);
                        report(spin_system,['  T[2, 0]:          sqrt(2/3)*(LzSz) x ' num2str(phi_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                        report(spin_system,['  T[2,-1]:                0.5*(LzSm) x ' num2str(phi_coupling(4)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    
                    % Prepare spherical tensor descriptors
                    nL(n,2,2)=L; nS(n,2,2)=S; opL(n,2,2)={'Lz'}; opS(n,2,2)={'Lz'};                                   ist_coeff{2}(n,2,2,3)=+sqrt(2/3);
                    nL(n,2,1)=L; nS(n,2,1)=S; opL(n,2,1)={'Lz'}; opS(n,2,1)={'L+'}; ist_coeff{1}(n,2,1,1)=+1/2;       ist_coeff{2}(n,2,1,2)=-1/2;
                    nL(n,2,3)=L; nS(n,2,3)=S; opL(n,2,3)={'Lz'}; opS(n,2,3)={'L-'}; ist_coeff{1}(n,2,3,3)=+1/2;       ist_coeff{2}(n,2,3,4)=+1/2;
                    
                case '*z'
                    
                    % Inform the user
                    report(spin_system,['(*,z) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    if norm(lam_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[1,+1]:               -0.5*(LpSz) x ' num2str(lam_coupling(1)/(2*pi),'%+0.5e') ' Hz']);
                        report(spin_system,['  T[1,-1]:               -0.5*(LmSz) x ' num2str(lam_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    if norm(phi_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[2,+1]:               -0.5*(LpSz) x ' num2str(phi_coupling(2)/(2*pi),'%+0.5e') ' Hz']);
                        report(spin_system,['  T[2, 0]:          sqrt(2/3)*(LzSz) x ' num2str(phi_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                        report(spin_system,['  T[2,-1]:                0.5*(LmSz) x ' num2str(phi_coupling(4)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    
                    % Prepare spherical tensor descriptors
                    nL(n,2,2)=L; nS(n,2,2)=S; opL(n,2,2)={'Lz'}; opS(n,2,2)={'Lz'};                                   ist_coeff{2}(n,2,2,3)=+sqrt(2/3);
                    nL(n,1,2)=L; nS(n,1,2)=S; opL(n,1,2)={'L+'}; opS(n,1,2)={'Lz'}; ist_coeff{1}(n,1,2,1)=-1/2;       ist_coeff{2}(n,1,2,2)=-1/2;
                    nL(n,3,2)=L; nS(n,3,2)=S; opL(n,3,2)={'L-'}; opS(n,3,2)={'Lz'}; ist_coeff{1}(n,3,2,3)=-1/2;       ist_coeff{2}(n,3,2,4)=+1/2;
                    
                case 'secular'
                    
                    % Inform the user
                    report(spin_system,['secular part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    if norm(lam_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[1, 0]:            -sqrt(1/8)*(LpSm-LmSp) x ' num2str(lam_coupling(2)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    if norm(phi_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[2, 0]: sqrt(2/3)*(LzSz-0.25*(LpSm+LmSp)) x ' num2str(phi_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    
                    % Prepare spherical tensor descriptors
                    nL(n,2,2)=L; nS(n,2,2)=S; opL(n,2,2)={'Lz'}; opS(n,2,2)={'Lz'};                                   ist_coeff{2}(n,2,2,3)=+sqrt(2/3);
                    nL(n,1,3)=L; nS(n,1,3)=S; opL(n,1,3)={'L+'}; opS(n,1,3)={'L-'}; ist_coeff{1}(n,1,3,2)=-sqrt(1/8); ist_coeff{2}(n,1,3,3)=-sqrt(2/3)/4;
                    nL(n,3,1)=L; nS(n,3,1)=S; opL(n,3,1)={'L-'}; opS(n,3,1)={'L+'}; ist_coeff{1}(n,3,1,2)=+sqrt(1/8); ist_coeff{2}(n,3,1,3)=-sqrt(2/3)/4;
                    
                case 'zz'
                    
                    % Inform the user
                    report(spin_system,['(z,z) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    if norm(phi_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[2, 0]:          sqrt(2/3)*(LzSz) x ' num2str(phi_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    
                    % Prepare spherical tensor descriptors
                    nL(n,2,2)=L; nS(n,2,2)=S; opL(n,2,2)={'Lz'}; opS(n,2,2)={'Lz'};                                   ist_coeff{2}(n,2,2,3)=+sqrt(2/3);
                    
                case 'z+'
                    
                    % Inform the user
                    report(spin_system,['(z,+) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    if norm(lam_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[1,+1]:          0.5*(LzSp) x ' num2str(lam_coupling(1)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    if norm(phi_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[2,+1]:         -0.5*(LzSp) x ' num2str(phi_coupling(2)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    
                    % Prepare spherical tensor descriptors
                    nL(n,2,1)=L; nS(n,2,1)=S; opL(n,2,1)={'Lz'}; opS(n,2,1)={'L+'}; ist_coeff{1}(n,2,1,1)=+1/2;       ist_coeff{2}(n,2,1,2)=-1/2;
                    
                case '+z'
                    
                    % Inform the user
                    report(spin_system,['(+,z) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    if norm(lam_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[1,+1]:         -0.5*(LpSz) x ' num2str(lam_coupling(1)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    if norm(phi_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[2,+1]:         -0.5*(LpSz) x ' num2str(phi_coupling(2)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    
                    % Prepare spherical tensor descriptors
                    nL(n,1,2)=L; nS(n,1,2)=S; opL(n,1,2)={'L+'}; opS(n,1,2)={'Lz'}; ist_coeff{1}(n,1,2,1)=-1/2;       ist_coeff{2}(n,1,2,2)=-1/2;
                    
                case 'T(L,+1)'
                    
                    % Inform the user
                    report(spin_system,['T(L,+1) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    if norm(lam_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[1,+1]:    -0.5*(LpSz-LzSp) x ' num2str(lam_coupling(1)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    if norm(phi_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[2,+1]:    -0.5*(LpSz+LzSp) x ' num2str(phi_coupling(2)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    
                    % Prepare spherical tensor descriptors
                    nL(n,2,1)=L; nS(n,2,1)=S; opL(n,2,1)={'Lz'}; opS(n,2,1)={'L+'}; ist_coeff{1}(n,2,1,1)=+1/2;       ist_coeff{2}(n,2,1,2)=-1/2;
                    nL(n,1,2)=L; nS(n,1,2)=S; opL(n,1,2)={'L+'}; opS(n,1,2)={'Lz'}; ist_coeff{1}(n,1,2,1)=-1/2;       ist_coeff{2}(n,1,2,2)=-1/2;
                    
                case 'z-'
                    
                    % Inform the user
                    report(spin_system,['(z,-) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    if norm(lam_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[1,-1]:          0.5*(LzSm) x ' num2str(lam_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    if norm(phi_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[2,-1]:          0.5*(LzSm) x ' num2str(phi_coupling(4)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    
                    % Prepare spherical tensor descriptors
                    nL(n,2,3)=L; nS(n,2,3)=S; opL(n,2,3)={'Lz'}; opS(n,2,3)={'L-'}; ist_coeff{1}(n,2,3,3)=+1/2;       ist_coeff{2}(n,2,3,4)=+1/2;
                    
                case '-z'
                    
                    % Inform the user
                    report(spin_system,['(-,z) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    if norm(lam_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[1,-1]:         -0.5*(LmSz) x ' num2str(lam_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    if norm(phi_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[2,-1]:          0.5*(LmSz) x ' num2str(phi_coupling(4)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    
                    % Prepare spherical tensor descriptors
                    nL(n,3,2)=L; nS(n,3,2)=S; opL(n,3,2)={'L-'}; opS(n,3,2)={'Lz'}; ist_coeff{1}(n,3,2,3)=-1/2;       ist_coeff{2}(n,3,2,4)=+1/2;
                    
                case 'T(L,-1)'
                    
                    % Inform the user
                    report(spin_system,['T(L,-1) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    if norm(lam_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[1,-1]:    -0.5*(LmSz-LzSm) x ' num2str(lam_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    if norm(phi_coupling,2)>spin_system.tols.liouv_zero
                        report(spin_system,['  T[2,-1]:     0.5*(LmSz+LzSm) x ' num2str(phi_coupling(4)/(2*pi),'%+0.5e') ' Hz']);
                    end
                    
                    % Prepare spherical tensor descriptors
                    nL(n,2,3)=L; nS(n,2,3)=S; opL(n,2,3)={'Lz'}; opS(n,2,3)={'L-'}; ist_coeff{1}(n,2,3,3)=+1/2;       ist_coeff{2}(n,2,3,4)=+1/2;
                    nL(n,3,2)=L; nS(n,3,2)=S; opL(n,3,2)={'L-'}; opS(n,3,2)={'Lz'}; ist_coeff{1}(n,3,2,3)=-1/2;       ist_coeff{2}(n,3,2,4)=+1/2;
                    
                case 'ignore'
                    
                    % Inform the user
                    report(spin_system,['anisotropic coupling ignored for spins ' num2str(L) ',' num2str(S) '.']);
                    
                otherwise
                    
                    % Bomb out with unexpected strength parameters
                    error(['unknown strength spec for the bilinear coupling between spins ' num2str(L) ' and ' num2str(S)]);
                    
            end
            
        end
        
    end
    
end
    
% Pack two-spin descriptor table
D2=table(reshape(nL, [9*size(pair_list,1) 1]),...
         reshape(nS, [9*size(pair_list,1) 1]),...
         reshape(opL,[9*size(pair_list,1) 1]),...
         reshape(opS,[9*size(pair_list,1) 1]),...
         reshape(isotropic,[9*size(pair_list,1) 1]),...
        [reshape(ist_coeff{1},[9*size(pair_list,1) 3])...
         reshape(ist_coeff{2},[9*size(pair_list,1) 5])],...
        [reshape(irr_comp{1}, [9*size(pair_list,1) 3])...
         reshape(irr_comp{2}, [9*size(pair_list,1) 5])],...
         'VariableNames',{'nL','nS','opL','opS','isotropic',...
                          'ist_coeff','irr_comp'});

% Kill insignificant rows
mask_a=(abs(D2.isotropic)       <spin_system.tols.liouv_zero);
mask_b=(sum(abs(D2.irr_comp),2) <spin_system.tols.liouv_zero);
mask_c=(sum(abs(D2.ist_coeff),2)<spin_system.tols.liouv_zero);
D2(mask_a&(mask_b|mask_c),:)=[];

% Merge descriptors 
descr=[D1; D2]; nterms=size(descr,1);
report(spin_system,[num2str(nterms) ' unique operators in the Hamiltonian descriptor.']);

% Tightest possible variable clean-up
clear('D1','D2','nL','nS','opL','opS','isotropic','ist_coeff',...
      'irr_comp','mask_a','mask_b','mask_c','L','S','coupling_iso',...
      'n','pair_list','quad_couplings');

% Balance the descriptor
descr=descr(randperm(size(descr,1)),:);

% Set up operator array
oper=cell(size(descr,1),1);

% Strip the spin system object down to minimum size
parfor_ss.sys=spin_system.sys; parfor_ss.tols=spin_system.tols;
parfor_ss.bas=spin_system.bas; parfor_ss.comp=spin_system.comp;

% Parfor rigging
if ~isworkernode
    D=parallel.pool.DataQueue;
    afterEach(D,@(~)parfor_progr);
    terms_done=0; last_toc=0;
    tic; ticBytes(gcp); do_diag=true;
else
    do_diag=false; D=[];
end

% Parfor progress updater
function parfor_progr()
    terms_done=terms_done+1; last_message=toc-last_toc;
    if (last_message>5)||(terms_done==nterms)
        report(spin_system,[num2str(terms_done) '/' num2str(nterms) ' operators done, ' ...
                            num2str(terms_done/toc) ' operators per second.']); 
        last_toc=toc;
    end
end

% Build component operators in XYZ form
report(spin_system,'building individual operators...'); tic;
parfor n=1:nterms
    
    % Get descriptor line
    descr_line=descr(n,:);
    
    % Compute operator from specification
    if descr_line.nS==0
        oper{n}=operator(parfor_ss,descr_line.opL,{descr_line.nL},operator_type,'xyz');
    else
        oper{n}=operator(parfor_ss,[descr_line.opL,descr_line.opS],...
                                   {descr_line.nL ,descr_line.nS },operator_type,'xyz');
    end

    % Report progress
    if do_diag, send(D,n); end

end

% Parfor communication stats
if ~isworkernode
    nbytes=mean(tocBytes(gcp),1)/2^20;
    report(spin_system,['average worker process received ' num2str(nbytes(1)) ...
                        ' MB and sent back ' num2str(nbytes(2)) ' MB']);
end

% Clean up and do the reporting
report(spin_system,['operator construction took ' num2str(toc()) ' seconds']);
report(spin_system,['operator array footprint: ' ...
                     num2str(whos('oper').bytes/1024^3) ' GB']);
clear('parfor_ss'); report(spin_system,'assembling the Hamiltonian...'); tic

% Get matrix dimension
dim=size(mprealloc(spin_system,0),1);

% Assemble isotropic part in XYZ form
H=cell(nterms,1);
for n=1:nterms
    if abs(descr.isotropic(n))>spin_system.tols.liouv_zero
        H{n}=[oper{n}(:,1) oper{n}(:,2) descr.isotropic(n)*oper{n}(:,3)];
        H{n}(abs(H{n}(:,3))<spin_system.tols.liouv_zero,:)=[];
    end
end
H(cellfun(@isempty,H))=[]; H=cell2mat(H);

% Convert XYZ form to sparse matrix
if isempty(H)
    H=mprealloc(spin_system,0);
else
    H=sparse(H(:,1),H(:,2),H(:,3),dim,dim); 
    H=clean_up(spin_system,H,spin_system.tols.liouv_zero);
end

% Assemble anisotropic part
if build_aniso

    % Rank loop
    for r=1:2

        % Preallocate output
        Q{r}=cell(2*r+1,2*r+1); %#ok<AGROW>

        % Projection loops
        for m=1:(2*r+1)
            for k=1:(2*r+1)

                % Assemble operator in XYZ form
                J=cell(nterms,1);
                for n=1:nterms

                    % Get the combined coefficient
                    coeff=descr.ist_coeff(n,r^2+k-1)*descr.irr_comp(n,r^2+m-1);

                    % Store the operator if significant
                    if abs(coeff)>spin_system.tols.liouv_zero
                        J{n}=[oper{n}(:,1) oper{n}(:,2) coeff*oper{n}(:,3)];
                        J{n}(abs(J{n}(:,3))<spin_system.tols.liouv_zero,:)=[];
                    end
                    
                end
                J(cellfun(@isempty,J))=[]; J=cell2mat(J);

                % Convert XYZ form to sparse matrix
                if isempty(J)
                    Q{r}{k,m}=mprealloc(spin_system,0);
                else
                    J=sparse(J(:,1),J(:,2),J(:,3),dim,dim);
                    Q{r}{k,m}=clean_up(spin_system,J,spin_system.tols.liouv_zero);
                end
                
            end
        end

    end

end

% Clean up and do the reporting
clear('J','descr','irr_comp','iso_coeff','oper');
report(spin_system,['Hamiltonian assembly took ' num2str(toc()) ' seconds']);

% Process giant spin Hamiltonian terms
if build_aniso

    % Detect the maximum spherical rank
    max_rank=max(cellfun(@numel,spin_system.inter.giant.coeff));
    report(spin_system,['maximum spherical rank in the giant spin Hamiltonian: ' num2str(max_rank)]);
    
    % Expand the data structure
    for r=1:max_rank
        if numel(Q)<r
            Q{r}=cell(2*r+1,2*r+1);
            Q{r}(1:(2*r+1),1:(2*r+1))={mprealloc(spin_system,0)};
        end
    end

    % Build the Hamiltonian
    for n=1:spin_system.comp.nspins
        
        % Loop over spherical ranks
        for r=1:numel(spin_system.inter.giant.coeff{n})
            
            switch spin_system.inter.giant.strength{n}
                
                case 'strong'
                   
                    % Inform the user
                    report(spin_system,['complete giant spin Hamiltonian of rank ' num2str(r) ' for spin ' num2str(n) '...']);
                    
                    % Loop over spherical tensor indices
                    for k=1:(2*r+1)
                        
                        % Build the tensor spec
                        ist_spec=['T' num2str(r) ',' num2str(r-k+1)];
                        
                        % Inform the user
                        report(spin_system,['           ' ist_spec ' x ' ...
                               num2str(spin_system.inter.giant.coeff{n}{r}(k)/(2*pi),'%+0.5e') ' Hz']);
                        
                        % Get the irreducible spherical tensor
                        T=operator(spin_system,{ist_spec},{n},operator_type);
                        
                        % Loop over coefficient indices
                        for m=1:(2*r+1)
                            
                            % Operator building
                            Q{r}{k,m}=Q{r}{k,m}+T*spin_system.inter.giant.coeff{n}{r}(m);
                            
                        end
                        
                    end
                    
                case 'secular'
                    
                    % Inform the user
                    report(spin_system,['secular giant spin Hamiltonian of rank ' num2str(r) ' for spin ' num2str(n) '...']);
                    
                    % Build the tensor spec
                    ist_spec=['T' num2str(r) ',' num2str(0)];
                    
                    % Inform the user
                    report(spin_system,['           ' ist_spec ' x ' ...
                           num2str(spin_system.inter.giant.coeff{n}{r}(r+1)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Only use the zero projection for T_lm
                    T=operator(spin_system,{ist_spec},{n},operator_type);
                    
                    % Loop over coefficient indices
                    for m=1:(2*r+1)
                        
                        % Operator building
                        Q{r}{r+1,m}=Q{r}{r+1,m}+T*spin_system.inter.giant.coeff{n}{r}(m);
                        
                    end
                    
                case 'ignore'
                    
                    % Tell the user
                    report(spin_system,['giant spin Hamiltonian terms of rank ' num2str(r)...
                                        ' for spin ' num2str(n) ' have been ignored.']); 
                    
            end
            
        end
                    
    end
    
end

% Remind the user about the anisotropic part
if ~build_aniso
    report(spin_system,'WARNING - only the isotropic part has been returned.');
end

% Warn about all-zero Hamiltonians
if (~build_aniso)&&(nnz(H)==0)
    report(spin_system,'WARNING - the Hamiltonian is ALL ZERO.');
end

% Print memory diagnostics
H_whos=whos('H'); report(spin_system,['memory footprint of H array: ' num2str(H_whos.bytes/1024^3) ' GB']);
if build_aniso
    Q_whos=whos('Q'); report(spin_system,['memory footprint of Q array: ' num2str(Q_whos.bytes/1024^3)  ' GB']);
end

end

% Consistency enforcement
function grumble(spin_system,operator_type)
if ~isfield(spin_system,'bas')
    error('basis set information is missing, run basis() before calling this function.');
end
if (~isfield(spin_system.inter.coupling,'strength'))||...
   (~isfield(spin_system.inter.zeeman,'strength'))     
    error('assumption information is missing, run assume() before calling this function.');
end
% Catch operator type errors
if (~ischar(operator_type))||...
   (~ismember(operator_type,{'left','right','comm','acomm'}))
    error('incorrect operator type specification.');
end
end

% IK's 2018 email to the European Research Council:
%
% Dear ERC,
% 
% I appreciate the panel invitation, but I must decline - making over a hundred
% skilled, persistent, brilliant and motivated enemies in the space of one week
% is not something I would be willing to do for any amount of money!  
% 
% It is also ironic that I should receive this invitation from an organisation
% that keeps rejecting my own proposals. This looks to me like being repeatedly 
% rejected for a date, and then asked by the same girl to share my views on the
% other guys she considers dating.
% 
% I am afraid not, and I would make the same suggestion as the one I had always
% been making in such cases since I was twelve: go ask somebody who you are ac-
% tually funding.
% 
% Yours sincerely,
% Ilya Kuprov.

