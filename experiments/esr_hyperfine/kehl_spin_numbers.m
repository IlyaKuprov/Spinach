% Spin quantum numbers from Spinach isotope labels. Syntax:
%
%      qnums=kehl_spin_numbers(labels)
%
% Parameters:
%
%   labels           - cell array of Spinach isotope labels.
%
% Outputs:
%
%   qnums            - spin quantum numbers.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_spin_numbers.m>

function qnums=kehl_spin_numbers(labels)

    % Check consistency
    grumble(labels);

    % Query Spinach isotope multiplicities
    qnums=zeros(numel(labels),1);
    for n=1:numel(labels)
        [~,multiplicity]=spin(labels{n});
        qnums(n)=(multiplicity-1)/2;
    end

end

% Consistency enforcement
function grumble(labels)
    if (~iscell(labels))||(~isvector(labels))
        error('labels must be a cell vector.');
    end
    if ~all(cellfun(@(x)ischar(x)||isstring(x),labels))
        error('all elements of labels must be character strings or string scalars.');
    end
end

