% Heuristic backend selector for step().
%
% The function returns a function handle from the supplied backend table. If
% the selected accelerated backend is unavailable or not beneficial, the
% default Spinach Taylor backend is returned.
%
% ilya.kuprov@weizmann.ac.il

function backend=step_heuristics(stats,backends)

% Cache path lookup for MATLAB's expmv()
persistent expmv_available

% Default to the native Spinach backend
backend=backends.default;

% Trivial cases should stay on the default path
if (stats.time_step==0)||(stats.dimension==0), return; end

% Use the scale already computed by step(); it is invariant under the
% Spinach sign convention A=-1i*L, and avoiding a second matrix norm keeps
% the selector cheap enough for short propagation calls.
L=stats.matrix;
alpha=stats.norm_mat;
if isa(alpha,'gpuArray'), alpha=gather(alpha); end
backend_name=local_backend_name(L,alpha);
if strcmp(backend_name,'vik'), return; end

% Return the selected backend handle
switch backend_name

    case 'expmv'

        % Fall back if Matlab expmv is unavailable
        if isempty(expmv_available)
            expmv_available=(exist('expmv','file')==2);
        end
        if expmv_available
            backend=backends.expmv;
        end

    case 'tay1'

        backend=backends.tay1;

    case 'tay2'

        backend=backends.tay2;

    otherwise

        backend=backends.default;

end

end

% Backend selector
function backend=local_backend_name(A,alpha)

% Get problem dimension
n=size(A,1);

% Dispatch by storage type, dimension, and scale
if isa(A,'gpuArray')

    % Dense GPU rule
    if alpha<=5
        backend='vik';
    elseif alpha<10
        if n>=1024
            backend='vik';
        else
            backend='tay1';
        end
    elseif alpha<50
        backend='tay2';
    elseif alpha<80
        backend='tay1';
    else
        if n<=256
            backend='tay2';
        else
            backend='tay1';
        end
    end

elseif issparse(A)

    % Sparse CPU rule
    if n>=4096
        if alpha<=25
            backend='expmv';
        elseif alpha<50
            backend='tay2';
        elseif alpha<80
            backend='tay1';
        else
            backend='tay2';
        end
    elseif n>=2048
        if alpha<15
            backend='expmv';
        elseif alpha<50
            backend='tay1';
        elseif alpha<80
            backend='tay2';
        else
            backend='tay1';
        end
    elseif n>=1024
        if alpha<4
            backend='vik';
        elseif alpha<10
            backend='expmv';
        elseif alpha<80
            backend='tay2';
        else
            backend='tay1';
        end
    elseif n>=512
        if alpha<5
            backend='vik';
        elseif alpha<9
            backend='expmv';
        else
            backend='tay2';
        end
    elseif n>=256
        if alpha<8
            backend='vik';
        elseif alpha<15
            backend='tay2';
        elseif alpha<80
            backend='tay1';
        else
            backend='tay2';
        end
    else
        if alpha<5
            backend='vik';
        elseif alpha<50
            backend='tay1';
        elseif alpha<80
            backend='tay2';
        else
            backend='tay1';
        end
    end

else

    % Dense CPU rule
    if n>=512
        if alpha<=8
            backend='vik';
        elseif alpha<50
            backend='expmv';
        else
            backend='tay1';
        end
    elseif n>=256
        if alpha<=6
            backend='vik';
        elseif alpha<10
            backend='expmv';
        elseif alpha<=12
            backend='vik';
        elseif alpha<25
            backend='expmv';
        elseif alpha<80
            backend='tay2';
        else
            backend='tay1';
        end
    elseif n>=128
        if alpha<=4
            backend='vik';
        elseif alpha<10
            backend='expmv';
        elseif alpha<15
            backend='tay2';
        elseif alpha<25
            backend='tay1';
        elseif alpha<50
            backend='tay2';
        else
            backend='tay1';
        end
    else
        if alpha<=5
            backend='vik';
        elseif alpha<50
            backend='tay1';
        else
            backend='tay2';
        end
    end

end

end
