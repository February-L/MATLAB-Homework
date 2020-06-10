function integral_val = int_simpson(integrand_func, lower_x, upper_x, step_len)
% Integrates a 1-dimension function (or more functions) by using the Simpson's method.
% Fu Chunhao, 4 June 2020
%
% function integral_val = int_simpson(integrand_func, lower_x, upper_x, step_len)
%     func_vals = integrand_func(lower_x: step_len: upper_x);
% function integral_val = int_simpson(integrand_func, x, step_len)
%     % For each i as the index of x, there should be (x(i + 1) - x(i) == step_len).
%     % x should be a row vector.
%     func_vals = integrand_func(x);
% function integral_val = int_simpson(integrand_func_vals, step_len)
%     % integrand_func_vals here should equal to the previous integrand_func(x).
    narginchk(2, 4);
    switch (nargin)
        case (2)
            integrand_func(:, [1, end]) = 0.5 .* integrand_func(:, [1, end]);
            integrand_func(:, 2: 2: (end - 1)) = 2.0 .* integrand_func(:, 2: 2: (end - 1));
            integral_val = 2.0/3.0 .* lower_x .* sum(integrand_func, 2);
        case (3)
            func_vals = integrand_func(lower_x);
            func_vals(:, [1, end]) = 0.5 .* func_vals(:, [1, end]);
            func_vals(:, 2: 2: (end - 1)) = 2.0 .* func_vals(:, 2: 2: (end - 1));
            integral_val = 2.0/3.0 .* upper_x .* sum(func_vals, 2);
        case (4)
            func_vals = integrand_func(lower_x: step_len: upper_x);
            func_vals(:, [1, end]) = 0.5 .* func_vals(:, [1, end]);
            func_vals(:, 2: 2: (end - 1)) = 2.0 .* func_vals(:, 2: 2: (end - 1));
            integral_val = 2.0/3.0 .* step_len .* sum(func_vals, 2);
    end
end
