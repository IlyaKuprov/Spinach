% Computes the transpose of a matrix in a polyadic representa-
% tion. Syntax:
%
%                         p=transpose(p)
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/transpose.m>

function p=transpose(p)

% Transpose every core
for n=1:numel(p.cores)
    for k=1:numel(p.cores{n})
        if isa(p.cores{n}{k},'function_handle')
            action=p.cores{n}{k}; adjoint=p.core_adj{n}{k};
            p.cores{n}{k}=@(x)conj(adjoint(conj(x)));
            p.core_adj{n}{k}=@(x)conj(action(conj(x)));
            p.core_dims{n}{k}=fliplr(p.core_dims{n}{k});
        else
            p.cores{n}{k}=transpose(p.cores{n}{k});
        end
    end
end

% Process prefix and suffix
new_prefix=fliplr(p.suffix);
for n=1:numel(new_prefix)
    new_prefix{n}=transpose(new_prefix{n});
end
new_suffix=fliplr(p.prefix);
for n=1:numel(new_suffix)
    new_suffix{n}=transpose(new_suffix{n});
end
p.prefix=new_prefix; p.suffix=new_suffix;

end


% Frantic orthodoxy is never rooted in faith 
% but in doubt. It is when we are unsure that
% we are doubly sure.
%
% Reinhold Niebuhr
% #NHEAD #NGRUM


