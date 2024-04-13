% Guesses a reasonable amide bond 15N CSA tensor anisotropy, and a 
% reasonable 13C=O tensor anisotropy, given a local protein geome-
% try. The tensors are oriented roughly according to
%
%         https://doi.org/10.1007/s10858-006-9037-6
%         https://doi.org/10.1021/ja00083a028
%         https://doi.org/10.1021/ja042863o
%
% Proline is not currently handled, amino acids are assumed to be 
% numbered from N-terminus to C-terminus.
%
% Note: these CSAs are very approximate. For accurate relaxation
%       analysis you must supply your own tensors.
%
% Note: this is an auxiliary function that is called by protein.m
%       protein import module. Direct calls are discouraged.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=Guess_csa_pro.m>

function CSAs=guess_csa_pro(aa_nums,pdb_ids,coords)

% Check consistency
grumble(aa_nums,pdb_ids,coords);

% Preallocate CSA array
CSAs=cell(numel(pdb_ids),1);

% Number the atoms
numbers=1:numel(coords);

% Loop over amino acids
for n=1:(max(aa_nums)-1)
    
    % Assign amide nitrogen CSAs
    if ismember('C',pdb_ids(aa_nums==n))&&...
       ismember('N',pdb_ids(aa_nums==(n+1)))&&...
       ismember('H',pdb_ids(aa_nums==(n+1)))
    
        % Get C coordinates
        local_coords=coords(aa_nums==n);
        C=local_coords{strcmp('C',pdb_ids(aa_nums==n))}; C=C(:);
   
        % Get N coordinates
        local_coords=coords(aa_nums==(n+1));
        N=local_coords{strcmp('N',pdb_ids(aa_nums==(n+1)))}; N=N(:);
   
        % Get H coordinates
        local_coords=coords(aa_nums==(n+1));
        H=local_coords{strcmp('H',pdb_ids(aa_nums==(n+1)))}; H=H(:);
        
        % Get the primary directions
        N_CO_vec=C-N; N_H_vec=H-N;
        
        % Double-check the distances
        if (norm(N_CO_vec,2)>2.0)||(norm(N_H_vec,2)>2.0)
            error('Amino acid numbering is not sequential.');
        end
        
        % Make ZZ eigenvector collinear with N-CO bond
        zz_eigvec=N_CO_vec;
        zz_eigvec=zz_eigvec/norm(zz_eigvec,2);
        
        % Make YY eigenvector perpendicular to the peptide plane
        yy_eigvec=cross(N_CO_vec,N_H_vec);
        yy_eigvec=yy_eigvec/norm(yy_eigvec,2);
        
        % Make XX eigenvector perpendicular to the other two
        xx_eigvec=cross(yy_eigvec,zz_eigvec);
        xx_eigvec=xx_eigvec/norm(xx_eigvec,2);
        
        % Build the eigenvalue matrix
        D=diag([-125 45 80]);
        
        % Build the eigenvector matrix
        V=[xx_eigvec yy_eigvec zz_eigvec];
        
        % Identify the nitrogen
        local_numbers=numbers(aa_nums==(n+1));
        nitrogen_number=local_numbers(strcmp('N',pdb_ids(aa_nums==(n+1))));
        
        % Compose the tensor
        CSAs{nitrogen_number}=V*D*V';
        
        % Report to the user
        disp(['Amide nitrogen CSA for residue ' num2str(n) ' guessed from local geometry.']);
    
    else

        % Report to the user
        disp(['CSA GUESS FAILED for the amide nitrogen in residue ' num2str(n)]);

    end

end

% Loop over amino acids
for n=1:(max(aa_nums)-1)
    
    % Assign C=O carbon CSAs
    if ismember('C',pdb_ids(aa_nums==n))&&...
       ismember('N',pdb_ids(aa_nums==(n+1)))&&...
       ismember('CA',pdb_ids(aa_nums==n))
    
        % Get C coordinates
        local_coords=coords(aa_nums==n);
        C=local_coords{strcmp('C',pdb_ids(aa_nums==n))}; C=C(:);
   
        % Get N coordinates
        local_coords=coords(aa_nums==(n+1));
        N=local_coords{strcmp('N',pdb_ids(aa_nums==(n+1)))}; N=N(:);
   
        % Get CA coordinates
        local_coords=coords(aa_nums==n);
        CA=local_coords{strcmp('CA',pdb_ids(aa_nums==n))}; CA=CA(:);
        
        % Get the primary directions
        N_C_vec=C-N; C_CA_vec=CA-C;
        
        % Double-check the distances
        if (norm(N_C_vec,2)>2.0)||(norm(C_CA_vec,2)>2.0)
            error('Amino acid numbering is not sequential.');
        end
        
        % Make XX eigenvector collinear with C-CA bond
        xx_eigvec=C_CA_vec;
        xx_eigvec=xx_eigvec/norm(xx_eigvec,2);
        
        % Make ZZ eigenvector perpendicular to the >C=O plane
        zz_eigvec=cross(C_CA_vec,N_C_vec);
        zz_eigvec=zz_eigvec/norm(zz_eigvec,2);
        
        % Make YY eigenvector perpendicular to the other two
        yy_eigvec=cross(xx_eigvec,zz_eigvec);
        yy_eigvec=yy_eigvec/norm(yy_eigvec,2);
        
        % Build the eigenvalue matrix
        D=diag([70  5  -75]);
        
        % Build the eigenvector matrix
        V=[xx_eigvec yy_eigvec zz_eigvec];
        
        % Identify the carbon
        local_numbers=numbers(aa_nums==n);
        carbon_number=local_numbers(strcmp('C',pdb_ids(aa_nums==n)));
        
        % Compose the tensor
        CSAs{carbon_number}=V*D*V';
        
        % Report to the user
        disp(['Carboxyl carbon CSA for residue ' num2str(n) ' guessed from local geometry.']);

    else

        % Report to the user
        disp(['CSA GUESS FAILED for the carboxyl carbon in residue ' num2str(n)]);

    end

end

% Loop over amino acids
for n=2:(max(aa_nums)-1)
    
    % Assign H-N proton CSAs
    if ismember('H',pdb_ids(aa_nums==n))&&...
       ismember('N',pdb_ids(aa_nums==n))&&...
       ismember('CA',pdb_ids(aa_nums==n))&&...
       ismember('C',pdb_ids(aa_nums==(n-1)))

        % Get C coordinates
        local_coords=coords(aa_nums==(n-1));
        C=local_coords{strcmp('C',pdb_ids(aa_nums==n))}; C=C(:);
    
        % Get H coordinates
        local_coords=coords(aa_nums==n);
        H=local_coords{strcmp('H',pdb_ids(aa_nums==n))}; H=H(:);
   
        % Get N coordinates
        local_coords=coords(aa_nums==n);
        N=local_coords{strcmp('N',pdb_ids(aa_nums==(n+1)))}; N=N(:);
   
        % Get CA coordinates
        local_coords=coords(aa_nums==n);
        CA=local_coords{strcmp('CA',pdb_ids(aa_nums==n))}; CA=CA(:);
        
        % Get the primary directions
        N_H_vec=N-H; N_CA_vec=N-CA; N_C_vec=N-C; H_CA_vec=H-CA;
        
        % Double-check the distances
        if (norm(N_H_vec,2)>1.2)||(norm(N_CA_vec,2)>1.6)||(norm(N_C_vec,2)>1.5)
            error('Amino acid numbering is not sequential.');
        end
        
        % Make YY eigenvector collinear with H-CA direction
        yy_eigvec=H_CA_vec;
        yy_eigvec=yy_eigvec/norm(yy_eigvec,2);
        
        % Make XX eigenvector perpendicular to the peptide plane
        xx_eigvec=cross(N_CA_vec,N_C_vec);
        xx_eigvec=xx_eigvec/norm(xx_eigvec,2);
        
        % Make ZZ eigenvector perpendicular to the other two
        zz_eigvec=cross(xx_eigvec,yy_eigvec);
        zz_eigvec=zz_eigvec/norm(zz_eigvec,2);
        
        % Build the eigenvalue matrix
        D=diag([7.0 0.0 -7.0]);
        
        % Build the eigenvector matrix
        V=[xx_eigvec yy_eigvec zz_eigvec];
        
        % Identify the proton
        local_numbers=numbers(aa_nums==n);
        proton_number=local_numbers(strcmp('H',pdb_ids(aa_nums==n)));
        
        % Compose the tensor
        CSAs{proton_number}=V*D*V';
        
        % Report to the user
        disp(['Amide proton CSA for residue ' num2str(n) ' guessed from local geometry.']);
    
    else

        % Report to the user
        disp(['CSA GUESS FAILED for the amide proton in residue ' num2str(n)]);

    end

end

end

% Consistency enforcement
function grumble(aa_nums,pdb_ids,coords)
if ~isnumeric(aa_nums)
    error('aa_nums must be a vector of integers.'); 
end
if ~iscell(pdb_ids)
    error('pdb_ids must be a cell array of character strings.');
end
if ~iscell(coords)
    error('coords must be a cell array of 3-vectors.');
end
if (numel(aa_nums)~=numel(pdb_ids))||(numel(pdb_ids)~=numel(coords))
    error('the input parameters must have the same number of elements.');
end
end

% "Where's the PCM solvent?"
%
% A 3rd year project student,
% rummaging through IK's sol-
% vent cabinet.

