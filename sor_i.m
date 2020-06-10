function x = sor_i(augmented_matrix, initial_guess, relaxation_factor, tolerance, max_step)
% Solves a system of linear equations by using the successive over-relaxation iteration (based on the Gauss-Seidel iteration).
% Note that not all invertible coefficient matrices lead to convergence.
% An appropriate value of the relaxation factor would be within interval (0, 2).
% Fu Chunhao, 24 April 2020
%
% function x = sor_i(augmented_matrix, initial_guess, relaxation_factor, tolerance, max_step)
%     x = initial_guess(:, 1);
%     x_count = length(x);
%     step_now = uint32(1); max_step = uint32(max_step);
%     while ((max_step == uint32(0)) || (step_now < max_step))
%         step_now = step_now + 1;
%         ... % (break by meeting the tolerance)
    x = initial_guess(:, 1);
    x_count = length(x);
    step_now = uint32(1); max_step = uint32(max_step);
    while ((max_step == uint32(0)) || (step_now < max_step))
        step_now = step_now + 1;
        tolerance_met = true;
        for i_row = 1: 1: x_count
            x_diff = (augmented_matrix(i_row, x_count + 1) - augmented_matrix(i_row, 1: 1: x_count) * x(:, 1)) ...
                     ./ augmented_matrix(i_row, i_row) .* relaxation_factor;
            x(i_row, 1) = x(i_row, 1) + x_diff;
            tolerance_met = tolerance_met & (abs(x_diff) < tolerance);
        end
        if (tolerance_met)
            break;
        end
    end
end
