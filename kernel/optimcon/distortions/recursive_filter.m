% Recursive filter also called Infinite Impulse Response.Applies digital filter based on 
% given coefficients to the input waveform.
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
%    a           - Feedback coefficients of the recursive filter, a vector [1, a1, a2]
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


function [y_filtered, J] = recursive_filter(w, a, b)

    % Initialize the filtered output with the same size as input
    y_filtered = zeros(size(w), 'like', w);

    % Number of channels and time slices
    [num_channels, num_cols] = size(w);

    % Initialize Jacobian as a sparse matrix
    if nargout > 1
        J = zeros(num_channels * num_cols, num_channels * num_cols);
    end

    % Iterate over each channel separately
    for c = 1:num_channels
        % Indices for the current channel in the Jacobian
        channel_start = (c-1)*num_cols + 1;
        channel_end = c*num_cols;

        % Initialize a temporary matrix to store derivatives for the current channel
        % This will be a lower triangular matrix
        if nargout > 1
            J_channel = zeros(num_cols, num_cols);
        end

        % Compute the filtered output and Jacobian for the current channel
        for n = 1:num_cols
            % Compute y_filtered(c, n) based on the recursive formula
            if n == 1
                y_filtered(c, n) = a(1) * w(c, n);
            elseif n == 2
                y_filtered(c, n) = a(1) * w(c, n) + a(2) * w(c, n-1) + b(1) * y_filtered(c, n-1);
            else
                y_filtered(c, n) = a(1) * w(c, n) + a(2) * w(c, n-1) + a(3) * w(c, n-2) ...
                                  + b(1) * y_filtered(c, n-1) + b(2) * y_filtered(c, n-2);
            end

            % Compute Jacobian entries if needed
            if nargout > 1

                % Initialize the current row in the Jacobian for y[n]
                J_current = zeros(1, num_cols);

                % Direct dependencies on inputs
                J_current(n) = a(1);  % d(y[n])/d(w[n]) = a0

                if n >= 2
                    J_current(n-1) = a(2);  % d(y[n])/d(w[n-1}) = a1
                end
                if n >= 3
                    J_current(n-2) = a(3);  % d(y[n])/d(w[n-2}) = a2
                end

                % Recursive dependencies on previous outputs
                if n >= 2
                    J_current = J_current + b(1) * J_channel(n-1, :);
                end
                if n >= 3
                    J_current = J_current + b(2) * J_channel(n-2, :);
                end

                % Assign the computed derivatives to the Jacobian matrix
                J_channel(n, :) = J_current;
            end
        end

        % Assign the channel-specific Jacobian to the overall Jacobian matrix
        if nargout > 1
            J(channel_start:channel_end, channel_start:channel_end) = J_channel;
        end
    end
end


