% Gather object from the GPU
%
% m.keitel@soton.ac.uk

function obj=gather(obj)

if obj.isGPU
    obj.row=gather(obj.row);
    obj.col=gather(obj.col);
    obj.val=gather(obj.val);
    obj.isGPU=false;
end

end

% Aerie, I've noticed the unfortunate fact that you live 
% by one of the great lessons of history that nothing is
% often a good thing to do and a clever thing to say.
%
% Edwin Odesseiron, in Baldur's Gate 2

