%KEHL_SPIN_NUMBERS Spin quantum numbers from Spinach isotope labels.
%
%   Spinach architecture migration May 2026 Talos

function qnums=kehl_spin_numbers(labels)
    qnums=zeros(numel(labels),1);
    for n=1:numel(labels)
        [~,multiplicity]=spin(labels{n});
        qnums(n)=(multiplicity-1)/2;
    end
end
