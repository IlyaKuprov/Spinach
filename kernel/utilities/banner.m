% Prints console banners. This is an internal function of 
% the kernel, user calls are discouraged. Syntax:
%
%              banner(spin_system,identifier)
%
% Parameters:
%
%    identifier - a character string with the banner name,
%                 see function text
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=banner.m>

function banner(spin_system,identifier)

% Check consistency
grumble(identifier);

% Print the banner
switch identifier
    case 'version_banner'
        report(spin_system,' ');
        report(spin_system,'============================================');
        report(spin_system,'=                                          =');
        report(spin_system,'=               SPINACH v2.11              =');
        report(spin_system,'=                                          =');
        report(spin_system,'=                <a href = "https://spindynamics.org/wiki/index.php?title=Spinach_developer_team">author list</a>               =');
        report(spin_system,'=                                          =');
        report(spin_system,'=               <a href = "https://spindynamics.org/wiki/index.php?title=Main_Page">documentation</a>              =');
        report(spin_system,'=                                          =');
        report(spin_system,'=                   <a href = "https://link.springer.com/book/10.1007/978-3-031-05607-9">book</a>                   =');
        report(spin_system,'=                                          =');
        report(spin_system,'=                MIT License               =');
        report(spin_system,'=                                          =');
        report(spin_system,'============================================');
        report(spin_system,' ');
    case 'spin_system_banner'
        report(spin_system,' ');
        report(spin_system,'============================================');
        report(spin_system,'=                                          =');
        report(spin_system,'=               SPIN SYSTEM                =');
        report(spin_system,'=                                          =');
        report(spin_system,'============================================');
        report(spin_system,' ');
    case 'basis_banner'
        report(spin_system,' ');
        report(spin_system,'============================================');
        report(spin_system,'=                                          =');
        report(spin_system,'=                BASIS SET                 =');
        report(spin_system,'=                                          =');
        report(spin_system,'============================================');
        report(spin_system,' ');
    case 'sequence_banner'
        report(spin_system,' ');
        report(spin_system,'============================================');
        report(spin_system,'=                                          =');
        report(spin_system,'=              PULSE SEQUENCE              =');
        report(spin_system,'=                                          =');
        report(spin_system,'============================================');
        report(spin_system,' ');
    case 'optimcon'
        report(spin_system,' ');
        report(spin_system,'============================================');
        report(spin_system,'=                                          =');
        report(spin_system,'=             OPTIMAL CONTROL              =');
        report(spin_system,'=                                          =');
        report(spin_system,'============================================');
        report(spin_system,' ');
    case 'optimisation'
        report(spin_system,' ');
        report(spin_system,'============================================');
        report(spin_system,'=                                          =');
        report(spin_system,'=               OPTIMISATION               =');
        report(spin_system,'=                                          =');
        report(spin_system,'============================================');
        report(spin_system,' ');
    otherwise
        error('unknown banner.');
end

end

% Consistency enforcement
function grumble(identifier)
if ~ischar(identifier)
    error('identifier must be a character string.');
end
end

% The free man will ask neither what his country can do 
% for him, nor what he can do for his country.
%
% Milton Friedman

