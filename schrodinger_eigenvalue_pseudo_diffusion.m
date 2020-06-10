clear;
% Calculating the eigen states of the stationary one-dimensional Schrodinger equation that the boundary values are given as 0:
potential_func = @(x) (0.0);
tolerance = 1e-8; sig_digits = int32(-log10(tolerance));
lower_x = 0.0; upper_x = 1.0;
x_step_len = 0.001;
x = lower_x: x_step_len: upper_x;
x_count = int32(length(x));
time_step_len = 0.0001;
to_the_time = 1.0;
least_time_step = int32(to_the_time / time_step_len);
    % This is to get eigen energy values that are more accurate than the
    % tolerance for highly non-degenerate states (the pseudo-diffusion
    % method converges fast for these states).
temp_factor = time_step_len / (x_step_len * x_step_len);
major_diag = zeros(1, x_count - 2) + 1.0 + 2.0 * temp_factor + time_step_len .* potential_func(x(2: 1: (x_count - 1)));
side_diag = zeros(1, x_count - 3) - temp_factor;
the_potential = potential_func(x);
global d2_dx2_the_wave; d2_dx2_the_wave = zeros(1, x_count);

to_nth_state = int32(5); % State 1 would be the ground state.
initial_guesses = {@(x) (x .* (1.0 - x)), ...
                   @(x) (x .* (x - 0.5) .* (1.0 - x)), ...
                   @(x) (x .* (x - 1.0/3.0) .* (x - 2.0/3.0) .* (1.0 - x)), ...
                   @(x) (x .* (x - 0.25) .* (x - 0.5) .* (x - 0.75) .* (1.0 - x)), ...
                   @(x) (x .* (x - 0.2) .* (x - 0.4) .* (x - 0.6) .* (x - 0.8) .* (1.0 - x))};
lower_states = zeros(to_nth_state, x_count);
fprintf("Step length of x: %G\nStep length of time: %G\nTheory        | Calculation   | Relative Error\n", x_step_len, time_step_len);
for i_state = 1: 1: to_nth_state
    state_factors = zeros(1, i_state - 1);
    last_wave = initial_guesses{i_state}(x);
    if (i_state >= int32(2))
        for i_state_inner = 1: 1: (i_state - 1)
            state_factors(i_state_inner) = int_trapezoidal(lower_states(i_state_inner, :) .* last_wave, x_step_len);
        end
        last_wave = last_wave - state_factors * lower_states(1: 1: (i_state - 1), :);
    end
    last_wave = last_wave ./ sqrt(int_trapezoidal(last_wave .* last_wave, x_step_len));
    calc_d2_dx2_the_wave(last_wave, x_step_len, x_count);
    eigen_energy = int_trapezoidal(last_wave .* (-d2_dx2_the_wave + the_potential .* last_wave), x_step_len);
    last_energy = eigen_energy - 10.0 * max(tolerance, eps(eigen_energy));
    the_wave = zeros(1, x_count) + 0.0;
    time_step_now = int32(0);
    while ((time_step_now < least_time_step) || (abs(eigen_energy - last_energy) > tolerance))
        time_step_now = time_step_now + 1;
        last_energy = eigen_energy;
        the_wave(2: 1: (x_count - 1)) = tridiag_chasing(side_diag, major_diag, side_diag, last_wave(2: 1: (x_count - 1)));
        if (i_state >= int32(2))
            for i_state_inner = 1: 1: (i_state - 1)
                state_factors(i_state_inner) = int_trapezoidal(lower_states(i_state_inner, :) .* the_wave, x_step_len);
            end
            the_wave = the_wave - state_factors * lower_states(1: 1: (i_state - 1), :);
        end
        the_wave = the_wave ./ sqrt(int_trapezoidal(the_wave .* the_wave, x_step_len));
        calc_d2_dx2_the_wave(the_wave, x_step_len, x_count);
        eigen_energy = int_trapezoidal(the_wave .* (-d2_dx2_the_wave + the_potential .* the_wave), x_step_len);
        last_wave = the_wave;
    end
    the_wave = the_wave ./ sqrt(int_simpson(the_wave .* the_wave, x_step_len));
    lower_states(i_state, :) = the_wave;
    calc_d2_dx2_the_wave(the_wave, x_step_len, x_count);
    eigen_energy = int_simpson(the_wave .* (-d2_dx2_the_wave + the_potential .* the_wave), x_step_len) ./ int_simpson(the_wave .* the_wave, x_step_len);
    actual_eigen_val = (round(sqrt(eigen_energy / (pi * pi))) * pi) ^ 2.0;
    fprintf("%#-13.*f | %#-13.*f | %#-.6f%%\n", sig_digits, actual_eigen_val, sig_digits, eigen_energy, 100.0 * (eigen_energy - actual_eigen_val) / actual_eigen_val);
    figure();
    plot(x, the_wave, "k-");
    xlabel("x"); ylabel("u(x)");
    switch (i_state)
        case (1)
            title("Î±À©É¢ - »ùÌ¬");
        case (2)
            title("Î±À©É¢ - µÚÒ»¼¤·¢Ì¬");
        case (3)
            title("Î±À©É¢ - µÚ¶þ¼¤·¢Ì¬");
        case (4)
            title("Î±À©É¢ - µÚÈý¼¤·¢Ì¬");
        case (5)
            title("Î±À©É¢ - µÚËÄ¼¤·¢Ì¬");
    end
end

function calc_d2_dx2_the_wave(the_wave, x_step_len, x_count)
    global d2_dx2_the_wave;
    for i_x = 2: 1: (x_count - 1)
        d2_dx2_the_wave(i_x) = (the_wave(i_x + 1) + the_wave(i_x - 1) - 2.0 * the_wave(i_x)) / (x_step_len * x_step_len);
    end
    d2_dx2_the_wave(1) = d2_dx2_the_wave(2);
    d2_dx2_the_wave(x_count) = d2_dx2_the_wave(x_count - 1);
end

% function integral_val = int_trapezoidal(integrand_func, lower_x, upper_x, step_len)
% function integral_val = int_simpson(integrand_func, lower_x, upper_x, step_len)
% function x = tridiag_chasing(lower_diagonal, major_diagonal, upper_diagonal, inhomogeneous)
