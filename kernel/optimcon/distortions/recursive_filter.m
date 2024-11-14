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


    % Number of channels and time slices
    [num_channels, num_cols] = size(w);

    % Initialize the filtered output
    y_filtered = zeros(size(w), 'like', w);

    % Initialize lists for Jacobian entries
    row_indices = [];
    col_indices = [];
    values = [];

    % Loop over channel pairs
    for k = 1:(num_channels/2)
        % Indices for current channel pair
        idx1 = 2*k - 1;
        idx2 = 2*k;

        % Extract u1 and u2
        u1 = w(idx1, :);
        u2 = w(idx2, :);

        % Compute amplitude and phase
        A = sqrt(u1.^2 + u2.^2);
        phi = atan2(u2, u1);

        % Handle zero amplitudes
        epsilon = 1e-12;
        A_safe = A;
        A_safe(A == 0) = epsilon;

        % Compute derivatives of A and phi
        dA_du1 = u1 ./ A_safe;
        dA_du2 = u2 ./ A_safe;
        dphi_du1 = -u2 ./ (A_safe.^2);
        dphi_du2 = u1 ./ (A_safe.^2);

        % Compute A_dist using filter
        b_coeffs = a(:).';
        a_coeffs = [1, -b(:).'];
        A_dist = filter(b_coeffs, a_coeffs, A);
        
        % Compute impulse response for sensitivity matrix
        impulse = zeros(1, num_cols);
        impulse(1) = 1;
        s_row = filter(b_coeffs, a_coeffs, impulse);
        
        % Ensure the first element of the row vector matches s_row(1)
        s_first_row = [s_row(1), zeros(1, num_cols - 1)];
        
        % Construct sensitivity matrix s without warning
        s = toeplitz(s_row, s_first_row);

        % Find all non-zero entries in s
        [n_indices, m_indices, s_values] = find(s);

        % Compute is_nm
        is_nm = (n_indices == m_indices);

        % Compute dy1_du1
        dy1_du1 = s_values.' .* dA_du1(m_indices) .* cos(phi(n_indices)) ...
                  - A_dist(n_indices) .* sin(phi(n_indices)) .* dphi_du1(m_indices) .* is_nm.';

        % Compute dy1_du2
        dy1_du2 = s_values.' .* dA_du2(m_indices) .* cos(phi(n_indices)) ...
                  - A_dist(n_indices) .* sin(phi(n_indices)) .* dphi_du2(m_indices) .* is_nm.';

        % Compute dy2_du1
        dy2_du1 = s_values.' .* dA_du1(m_indices) .* sin(phi(n_indices)) ...
                  + A_dist(n_indices) .* cos(phi(n_indices)) .* dphi_du1(m_indices) .* is_nm.';

        % Compute dy2_du2
        dy2_du2 = s_values.' .* dA_du2(m_indices) .* sin(phi(n_indices)) ...
                  + A_dist(n_indices) .* cos(phi(n_indices)) .* dphi_du2(m_indices) .* is_nm.';

        % Compute Jacobian indices
        y_idx1 = (n_indices - 1) * num_channels + idx1;
        y_idx2 = (n_indices - 1) * num_channels + idx2;
        w_idx1 = (m_indices - 1) * num_channels + idx1;
        w_idx2 = (m_indices - 1) * num_channels + idx2;

        % Accumulate indices and values
        row_indices = [row_indices; y_idx1; y_idx1; y_idx2; y_idx2];
        col_indices = [col_indices; w_idx1; w_idx2; w_idx1; w_idx2];
        values = [values; dy1_du1.'; dy1_du2.'; dy2_du1.'; dy2_du2.'];

        % Compute filtered output
        y_filtered(idx1, :) = A_dist .* cos(phi);
        y_filtered(idx2, :) = A_dist .* sin(phi);
    end

    % Construct the sparse Jacobian matrix
    Jacobian = sparse(row_indices, col_indices, values, numel(y_filtered), numel(w));

end


