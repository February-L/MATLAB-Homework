clear;
% Calculating the eigen states of the stationary one-dimensional Schrodinger equation that the boundary values are given as 0:
potential_func = @(x) (0.0);
tolerance = 1e-8; sig_digits = int32(-log10(tolerance));
lower_x = 0.0; upper_x = 1.0;
first_guess = 0.0;
guess_step_len = 5.0;
to_nth_state = int32(5);

for x_step_len = [0.02, 0.01, 0.005, 0.001]
    x = lower_x: x_step_len: upper_x;
    x_count = int32(length(x));
    appeared_eigen_vals = zeros(1, to_nth_state);
    eigen_count_now = int32(0);
    fprintf("Step length of x: %G. Using ", x_step_len);
    if (x_step_len >= 0.02)
        % The results of symmtridiag_eigen and inv_power_i are actually the same.
        % The accuracy of the calculated eigen values only depends on the
        % step length of x as long as the tolerance is small enough.
        fprintf("symmtridiag_eigen");
        side_diag = zeros(1, x_count - 1) - 1.0 / (x_step_len * x_step_len);
        major_diag = zeros(1, x_count) + 2.0 / (x_step_len * x_step_len) + potential_func(x);
    else
        fprintf("inv_power_i");
        the_matrix = diag((2.0 / (x_step_len * x_step_len) + potential_func(x)) .* ones(1, x_count));
        for i_side = 1: 1: (x_count - 1)
            the_matrix(i_side, i_side + 1) = -1.0 / (x_step_len * x_step_len);
            the_matrix(i_side + 1, i_side) = -1.0 / (x_step_len * x_step_len);
        end
        guess_vector = ones(x_count, 1);
    end
    fprintf(":\nNearest Theory |   Calculation | Relative Error\n");
    the_guess = first_guess;
    while (eigen_count_now < to_nth_state)
        the_guess = the_guess + guess_step_len;
        if (x_step_len >= 0.02)
            eigen_val = symmtridiag_eigen(major_diag, side_diag, the_guess, 0.5 * guess_step_len, tolerance, 0);
        else
            [eigen_vector, eigen_val] = inv_power_i(the_matrix, guess_vector, the_guess, tolerance, 50);
        end
        if (all(abs(eigen_val - appeared_eigen_vals) > abs(10.0 * x_step_len * eigen_val)))
            eigen_count_now = eigen_count_now + 1;
            appeared_eigen_vals(eigen_count_now) = eigen_val;
            actual_eigen_val = (round(sqrt(eigen_val / (pi * pi))) * pi) ^ 2.0;
            fprintf("%# 14.*f | %# 13.*f | %# 13.4f%%\n", sig_digits, actual_eigen_val, sig_digits, eigen_val, 100.0 * (eigen_val - actual_eigen_val) / actual_eigen_val);
            if (x_step_len == 0.001)
                eigen_vector = eigen_vector ./ sqrt(int_simpson((eigen_vector .* eigen_vector)', x_step_len));
                figure();
                plot(x, eigen_vector, "k-");
                xlabel("x"); ylabel("u(x)");
                switch (eigen_count_now)
                    case int32(1)
                        title("基态");
                    case int32(2)
                        title("第一激发态");
                    case int32(3)
                        title("第二激发态");
                    case int32(4)
                        title("第三激发态");
                    case int32(5)
                        title("第四激发态");
                end
            end
        end
    end
    fprintf("\n");
end

% 结果的准确性与初始猜测无关，故注释掉下面的部分。
% fprintf("\nStep length of x: %G. Using inv_power_i with a good guess of the ground state:\n", x_step_len);
% good_guess_vector = (x .* (1.0 - x))';
% [eigen_vector, eigen_val] = inv_power_i(the_matrix, good_guess_vector, 0.0, tolerance, 0);
% actual_eigen_val = (round(sqrt(eigen_val / (pi * pi))) * pi) ^ 2.0;
% fprintf("Calculated Energy: %# 15.*f\nTheoretical Energy: %# 14.*f\n", sig_digits + 1, eigen_val, sig_digits + 1, actual_eigen_val);
% eigen_vector = eigen_vector ./ sqrt(int_simpson((eigen_vector .* eigen_vector)', x_step_len));
% figure();
% plot(x, eigen_vector, "k-");
% xlabel("x"); ylabel("u(x)");
% title("基态（使用较好的初始猜测）");

% function integral_val = int_simpson(integrand_func, lower_x, upper_x, step_len)
% function eigen_val = symmtridiag_eigen(major_diagonal, side_diagonal, guess_eigen_val, init_step_len, tolerance, max_step)
% function [eigenvector, eigenvalue] = inv_power_i(the_matrix, guess_eigenvector, guess_eigenvalue, eigenvalue_tolerance, max_step)
