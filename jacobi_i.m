function x = jacobi_i(augmented_matrix, initial_guess, tolerance, max_step)
% Solves a system of linear equations by using the Jacobi iteration.
% Note that not all invertible coefficient matrices lead to convergence.
% Theorem: If the coefficient matrix is strictly diagonally dominant,
%          the Jacobian iteration converges on a unique solution.
% Fu Chunhao, 24 April 2020
%
% function x = jacobi_i(augmented_matrix, initial_guess, tolerance, max_step)
%     x = initial_guess(:, 1);
%     x_count = length(x);
%     x_diff = zeros(x_count, 1);
%     x_factor = diag(augmented_matrix);
%     step_now = uint32(1); max_step = uint32(max_step);
%     while ((max_step == uint32(0)) || (step_now < max_step))
%         step_now = step_now + 1;
%         ... % (break by meeting the tolerance)
    x = initial_guess(:, 1);
    x_count = length(x);
    x_diff = zeros(x_count, 1);
    x_factor = diag(augmented_matrix);
    step_now = uint32(1); max_step = uint32(max_step);
    while ((max_step == uint32(0)) || (step_now < max_step))
        step_now = step_now + 1;
        x_diff(:, 1) = (augmented_matrix(1: 1: x_count, x_count + 1) - augmented_matrix(1: 1: x_count, 1: 1: x_count) * x(:, 1)) ...
                       ./ x_factor(:, 1);
        x = x + x_diff;
        if (all(abs(x_diff) < tolerance))
            break;
        end
    end
end
