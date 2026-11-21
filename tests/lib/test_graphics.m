% Isolates regression-test graphics until the returned guard is cleared.
% Syntax:
%
%                       cleaner=test_graphics()
%
% Outputs:
%
%     cleaner - onCleanup guard that deletes new figures and restores
%               root defaults and the previously current figure
%
% Existing figures are hidden from implicit plotting by clearing the root
% current figure. Tests must not explicitly delete caller-owned figures.

function cleaner=test_graphics()

% Save caller graphics state before imposing invisible test figures
figures=findall(groot,'Type','figure');
current=get(groot,'CurrentFigure');
defaults=get(groot,'default');
cleaner=onCleanup(@()restore_graphics(figures,current,defaults));
set(groot,'DefaultFigureVisible','off','CurrentFigure',[]);

end

% Restore graphics on normal return and during exception unwinding
function restore_graphics(figures,current,defaults)

% Delete only figures created inside this guard, including hidden handles
created=findall(groot,'Type','figure');
delete(created(~ismember(created,figures)));

% Remove test-added defaults before reinstating the caller's defaults
changed=fieldnames(get(groot,'default'));
for n=1:numel(changed)
    set(groot,changed{n},'remove');
end
set(groot,defaults);
if isempty(current)||isgraphics(current)
    set(groot,'CurrentFigure',current);
end

end


