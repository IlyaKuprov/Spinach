% Docs header here.

function test=alpha_conds(test_type,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system)

% Select test type
if test_type==0
    
    % Gradient free test - ensure increase
    test=(fx_1 > fx_0);
    
elseif test_type==1
    
    % Armijo condition - ensure sufficient increase
    test=(fx_1 >= fx_0 + spin_system.control.ls_c1*alpha*(gfx_0'*dir));
    
elseif test_type==2
    
    % Strong Wolfe (curvature) condition - avoid underestimated step length
    test=(abs(gfx_1'*dir) <= spin_system.control.ls_c2*abs(gfx_0'*dir));
    
elseif test_type==3
    
    % Gradient test - ensure an ascent direction
    test=(gfx_1'*dir > 0);
    
end

end

% These are illusions of popular history which 
% a successful religion must promote: evil men
% never prosper; only the brave deserve the fair;
% honesty is the best policy; actions speak lou-
% der than words; virtue always triumphs; a good 
% deed is its own reward; any bad human can be 
% reformed; religious talismans protect one from
% demon possession; only females understand the
% ancient mysteries; the rich are doomed to un-
% happiness...
%
% Frank Herbert, in the Dune series

