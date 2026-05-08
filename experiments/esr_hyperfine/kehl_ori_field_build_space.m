%KEHL_ORI_FIELD_BUILD_SPACE Magnetic-quantum-number basis table for EPR nuclei.
%
%   Spinach architecture migration May 2026 Talos

function M=kehl_ori_field_build_space(S)

% Check consistency
grumble(S);

    % builds up the space of basis vectors in columns,
    % script by M. Bennati
    %
    % input parameters:
    % S: spin quantum number
    %
    % output parameters:
    % M: basis vectors in columns
    %
    % February 2024 A. Kehl (akehl@gwdg.de)
    %

    % Initialize the space to a null
    M=[];

    % Initialize the dimension to one
    dim=1;

    % Loop over each spin
    for j=1:length(S)

      % Initialize the holder
      temp=[];

      % Loop over the z components
      for k=-S(j):S(j)

        %   adding block of the old space
        temp=[temp; M k*ones(dim,1)];

      %   plus a column of the z component
      end

      % Update M
      M=temp;

      % The dimensionality is the number
      dim=size(M,1);
                                            %   of rows in M
    end
end

% Consistency enforcement
function grumble(S)
if ~isnumeric(S)
    error('S must be numeric.');
end
end
