function B0G=coil_map_lookup(B0,BXMap,BYMap,BZMap,currents,xyz)

% Build interpolants in Angstroms, assume Tesla/Ampere for fields
F_B1_11=scatteredInterpolant(1e10*BXMap(:,1),1e10*BXMap(:,2),1e10*BXMap(:,3),...
                             BXMap(:,4),'linear','none');
F_B1_12=scatteredInterpolant(1e10*BXMap(:,1),1e10*BXMap(:,2),1e10*BXMap(:,3),...
                             BXMap(:,5),'linear','none');
F_B1_13=scatteredInterpolant(1e10*BXMap(:,1),1e10*BXMap(:,2),1e10*BXMap(:,3),...
                             BXMap(:,6),'linear','none');
F_B1_21=scatteredInterpolant(1e10*BYMap(:,1),1e10*BYMap(:,2),1e10*BYMap(:,3),...
                             BYMap(:,4),'linear','none');
F_B1_22=scatteredInterpolant(1e10*BYMap(:,1),1e10*BYMap(:,2),1e10*BYMap(:,3),...
                             BYMap(:,5),'linear','none');
F_B1_23=scatteredInterpolant(1e10*BYMap(:,1),1e10*BYMap(:,2),1e10*BYMap(:,3),...
                             BYMap(:,6),'linear','none');
F_B1_31=scatteredInterpolant(1e10*BZMap(:,1),1e10*BZMap(:,2),1e10*BZMap(:,3),...
                             BZMap(:,4),'linear','none');
F_B1_32=scatteredInterpolant(1e10*BZMap(:,1),1e10*BZMap(:,2),1e10*BZMap(:,3),...
                             BZMap(:,5),'linear','none');
F_B1_33=scatteredInterpolant(1e10*BZMap(:,1),1e10*BZMap(:,2),1e10*BZMap(:,3),...
                             BZMap(:,6),'linear','none');

% Preallocate output
B0G=cell(size(xyz)); B0G(:)={B0};

% Loop over coordinates
for n=1:numel(xyz)

    % X coil fields
    B0G{n}(1)=B0G{n}(1)+currents(1)*F_B1_11(xyz{n}(1),xyz{n}(2),xyz{n}(3)); 
    B0G{n}(2)=B0G{n}(2)+currents(1)*F_B1_12(xyz{n}(1),xyz{n}(2),xyz{n}(3));
    B0G{n}(3)=B0G{n}(3)+currents(1)*F_B1_13(xyz{n}(1),xyz{n}(2),xyz{n}(3));

    % Y coil fields
    B0G{n}(1)=B0G{n}(1)+currents(2)*F_B1_21(xyz{n}(1),xyz{n}(2),xyz{n}(3));
    B0G{n}(2)=B0G{n}(2)+currents(2)*F_B1_22(xyz{n}(1),xyz{n}(2),xyz{n}(3));
    B0G{n}(3)=B0G{n}(3)+currents(2)*F_B1_23(xyz{n}(1),xyz{n}(2),xyz{n}(3));

    % Z coil fields
    B0G{n}(1)=B0G{n}(1)+currents(3)*F_B1_31(xyz{n}(1),xyz{n}(2),xyz{n}(3));
    B0G{n}(2)=B0G{n}(2)+currents(3)*F_B1_32(xyz{n}(1),xyz{n}(2),xyz{n}(3));
    B0G{n}(3)=B0G{n}(3)+currents(3)*F_B1_33(xyz{n}(1),xyz{n}(2),xyz{n}(3));

end

end