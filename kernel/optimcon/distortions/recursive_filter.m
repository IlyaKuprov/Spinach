% Recursive filter also called Infinite Impulse Response. Applies digital filter based on 
% given coefficients to the amplitude of the input waveform in XY channels.
%
% The recursive filter formula is:
%
%                 y[n] = a1*w[n] +a2 * w[n-1] + a3 * w[n-2] + b1 * w[n-1] + b2 * w[n-2]
%
% Outputs the filtered waveform in the same layout as the input and the
% Jacobian matrix for the input array.
%
% Syntax:
%
%                 [y_filtered, J] = recursive_filter(w, a, b)
%
% Parameters:
%
%    w           - Input waveform, one time slice per column, rows arranged as
%                  separate channels (e.g., Channel1_X, Channel1_Y, Channel2_X, ...)
%
%    a           - Feedback coefficients of the recursive filter, a vector [a1, a2, a3]
%    b           - Feedforward coefficients of the recursive filter, a vector [b1,b2]
%
% Outputs:
%
%    y_filtered  - Filtered waveform in the same layout as the input
%    J           - Filter Jacobian matrix with respect to the input, sparse
%
% Example:
%
%   Single Pole Low Pass Filter:
%
%    a0=1-x
%    b1=x
%
%    a = [0.15, 0, 0];  % Example feedback coefficients
%    b = [0.85, 0];  % Example feedforward coefficients
%
%   Single Pole High Pass Filter:
%    
%    a0= 1+x / 2
%    a1= -(1+x) / 2
%    b1=x
%
%    a = [0.93, -0.93, 0];  % Example feedback coefficients
%    b = [0.86, 0];  % Example feedforward coefficients
%
%    Physically, x is the amount of decay between
%    adjacent samples
%
%    One is able to apply any general filter including a narrow-band filter
%    representing RLC response. However, this will depend on case by case
%    basis because coefficients have to be recalculated.
%
%    The Scientist and Engineer's Guide to
%    Digital Signal Processing (2009) Chapter 19
%    Second Edition
%    by  Steven W. Smith
%
%    https://ia801301.us.archive.org/23/items/GuideToDigitalSignalProcessing/Guide%20To%20Digital%20Signal%20Processing.pdf
%
%    [filtered_y, Jacobian] = recursive_filter(input_x, a, b);
%
% Author:
%    u.rasulov@soton.ac.uk
%


function [y_filtered, Jacobian] = recursive_filter(w, a, b)

    % Initialize the filtered output with the same size as input
    y_filtered = zeros(size(w), 'like', w);
    % Initialize the Jacobian matrix
    Jacobian = zeros(numel(y_filtered), numel(w), 'like', w);

    % Number of channels and time slices
    [num_channels, num_cols] = size(w);

    % Loop over channels
    for k = 1:(num_channels/2)
        % Indices for current channel
        idx1 = 2*k - 1;
        idx2 = 2*k;

        % Get amplitude and phase
        u1 = w(idx1, :);
        u2 = w(idx2, :);
        A = sqrt(u1.^2 + u2.^2);
        phi = atan2(u2, u1);

        % Initialize A_dist and its Jacobian
        A_dist = zeros(size(A), 'like', A);
        s = zeros(num_cols, num_cols, 'like', A);  % For storing derivatives

        % Compute A_dist and its Jacobian recursively
        for n = 1:num_cols
            % Compute A_dist(n)
            if n == 1
                A_dist(n) = a(1) * A(n);
            elseif n == 2
                A_dist(n) = a(1) * A(n) + a(2) * A(n-1) + b(1) * A_dist(n-1);
            else
                A_dist(n) = a(1) * A(n) + a(2) * A(n-1) + a(3) * A(n-2) ...
                            + b(1) * A_dist(n-1) + b(2) * A_dist(n-2);
            end
            % Compute s(n, m)
            for m = 1:n
                % Initialize s(n, m)
                s(n, m) = 0;
                % Direct dependencies
                if n == m
                    s(n, m) = s(n, m) + a(1);
                end
                if n - 1 == m
                    s(n, m) = s(n, m) + a(2);
                end
                if n - 2 == m
                    s(n, m) = s(n, m) + a(3);
                end
                % Recursive contributions
                if n > 1
                    s(n, m) = s(n, m) + b(1) * s(n - 1, m);
                end
                if n > 2
                    s(n, m) = s(n, m) + b(2) * s(n - 2, m);
                end
            end
        end

        % Compute derivatives of A and phi with respect to u1 and u2
        dA_du1 = u1 ./ A;
        dA_du2 = u2 ./ A;
        dphi_du1 = -u2 ./ (A.^2);
        dphi_du2 = u1 ./ (A.^2);

        % Compute the Jacobian entries
        for n = 1:num_cols
            for m = 1:n
                if s(n, m) ~= 0
                    % Partial derivatives with respect to u1(m)
                    dy1_du1 = s(n, m) * dA_du1(m) * cos(phi(n));
                    dy1_du1 = dy1_du1 - A_dist(n) * sin(phi(n)) * dphi_du1(m) * (n == m);
                    dy2_du1 = s(n, m) * dA_du1(m) * sin(phi(n));
                    dy2_du1 = dy2_du1 + A_dist(n) * cos(phi(n)) * dphi_du1(m) * (n == m);

                    % Partial derivatives with respect to u2(m)
                    dy1_du2 = s(n, m) * dA_du2(m) * cos(phi(n));
                    dy1_du2 = dy1_du2 - A_dist(n) * sin(phi(n)) * dphi_du2(m) * (n == m);
                    dy2_du2 = s(n, m) * dA_du2(m) * sin(phi(n));
                    dy2_du2 = dy2_du2 + A_dist(n) * cos(phi(n)) * dphi_du2(m) * (n == m);

                    % Adjusted Indices for desired Jacobian structure
                    % y_idx corresponds to y_filtered(i,n)
                    % w_idx corresponds to w(k,m)
                    y_idx1 = (n - 1) * num_channels + idx1;
                    y_idx2 = (n - 1) * num_channels + idx2;
                    w_idx1 = (m - 1) * num_channels + idx1;
                    w_idx2 = (m - 1) * num_channels + idx2;

                    % Assign to Jacobian matrix
                    Jacobian(y_idx1, w_idx1) = Jacobian(y_idx1, w_idx1) + dy1_du1;
                    Jacobian(y_idx1, w_idx2) = Jacobian(y_idx1, w_idx2) + dy1_du2;
                    Jacobian(y_idx2, w_idx1) = Jacobian(y_idx2, w_idx1) + dy2_du1;
                    Jacobian(y_idx2, w_idx2) = Jacobian(y_idx2, w_idx2) + dy2_du2;
                end
            end
        end

        % Compute filtered output
        y_filtered(idx1, :) = A_dist .* cos(phi);
        y_filtered(idx2, :) = A_dist .* sin(phi);
    end
end


