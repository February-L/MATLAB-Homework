clear;
boundary_val = 0.0;
% Calculating the eigen states of the stationary one-dimensional Schrodinger equation that the boundary values are given as 0:
use_nonzero_potential_eg = false;
potential_func = @(x) (0.0);
if (use_nonzero_potential_eg)
    potential_func = @(x) (100.0.*exp(-10.0.*(x-0.5).^2.0));
end
e_tolerance = 1e-8; sig_digits = int32(-log10(e_tolerance));
lower_x = 0.0; upper_x = 1.0;
x_step_len = 0.001;
side_temp = x_step_len * x_step_len / 12.0;
centre_temp = x_step_len * x_step_len * 5.0 / 6.0;
x = lower_x: x_step_len: upper_x;
x_count = int32(length(x));
the_potential = zeros(1, x_count) + potential_func(x); 
the_wave = zeros(1, x_count) + boundary_val;

to_nth_state = int32(5); % State 1 would be the ground state.
i_state = int32(0);
eigen_energy = 0.0;
usual_energy_step_length = 1.0;
if (use_nonzero_potential_eg)
    fprintf("Step length of x: %G\n n |  Eigen Energy\n", x_step_len);
else
    fprintf("Step length of x: %G\nNearest Theory |   Calculation | Relative Error\n", x_step_len);
end
deri_diff_tolerance = 10.0 * x_step_len * x_step_len; % The 'deri' here means 'derivative'.
zero_tolerance = deri_diff_tolerance * x_step_len;
while (i_state < to_nth_state)
    i_state = i_state + 1;
    e_step_len = usual_energy_step_length;
    not_first_try = false;
    first_diff = 0.0;
    the_diff = 0.0;
    was_comparing_deri = false; % The 'deri' here means 'derivative'.
    while (true)
        x_trans_i = x_count;
        potential_e_diff = the_potential - eigen_energy;
        for i_x = 4: 1: (x_count - 2)
            if ((potential_e_diff(i_x) * potential_e_diff(i_x - 1) <= 0.0) && ((potential_e_diff(i_x) ~= 0.0) || (potential_e_diff(i_x - 1) ~= 0.0)))
                if (abs(potential_e_diff(i_x) - eigen_energy) < abs(potential_e_diff(i_x - 1) - eigen_energy))
                    x_trans_i = i_x;
                else
                    x_trans_i = i_x - 1;
                end
                break;
            end
        end
        if (x_trans_i ~= x_count)
            if ((~was_comparing_deri) && not_first_try)
                not_first_try = false;
                eigen_energy = eigen_energy - e_step_len;
            end
        else
            if (was_comparing_deri && not_first_try)
                not_first_try = false;
                eigen_energy = eigen_energy - e_step_len;
                x_trans_i = x_count;
            end
        end
        forward_wave_end = 0.0;
        for motion_of_i = [int32(1), int32(-1)]
            if (motion_of_i == int32(1))
                i_x = int32(2);
            else
                i_x = x_count - 1;
            end
            the_wave(i_x) = boundary_val + x_step_len * (1.0 + sqrt(abs(eigen_energy)));
            while (i_x ~= x_trans_i)
                the_wave(i_x + motion_of_i) = ( ...
                                                  (2.0 * the_wave(i_x) - the_wave(i_x - motion_of_i)) + ...
                                                  ( ...
                                                      (centre_temp * the_potential(i_x) * the_wave(i_x) + side_temp * the_potential(i_x - motion_of_i) * the_wave(i_x - motion_of_i)) - ...
                                                      (centre_temp * eigen_energy * the_wave(i_x) + side_temp * eigen_energy * the_wave(i_x - motion_of_i)) ...
                                                  ) ...
                                              ) / (1.0 - side_temp * (the_potential(i_x + motion_of_i) - eigen_energy));
                i_x = i_x + motion_of_i;
            end
            if (motion_of_i == int32(1))
                forward_wave_end = the_wave(x_trans_i);
                if (x_trans_i == x_count)
                    the_wave(x_count) = boundary_val;
                    break;
                end
            end
        end
        if (x_trans_i ~= x_count)
            deri_diff_tolerance = abs( ...
                                      0.5 * x_step_len * potential_e_diff(x_trans_i) * forward_wave_end + 1.0/6.0 * x_step_len * ( ...
                                          0.5 * (the_potential(x_trans_i + 1) - the_potential(x_trans_i - 1)) * forward_wave_end + ...
                                          potential_e_diff(x_trans_i) * (forward_wave_end - the_wave(x_trans_i - 1)) ...
                                      ) ...
                                  ) + x_step_len * x_step_len * x_step_len;
            zero_tolerance = deri_diff_tolerance * x_step_len;
            if (abs(forward_wave_end) > zero_tolerance)
                the_wave(1: 1: (x_trans_i - 1)) = the_wave(1: 1: (x_trans_i - 1)) .* (the_wave(x_trans_i) / forward_wave_end);
                forward_wave_end = the_wave(x_trans_i);
            else
                forward_wave_extension = ( ...
                                             (2.0 * forward_wave_end - the_wave(x_trans_i - 1)) + ...
                                             ( ...
                                                 (centre_temp * the_potential(x_trans_i) * forward_wave_end + side_temp * the_potential(x_trans_i - 1) * the_wave(x_trans_i - 1)) - ...
                                                 (centre_temp * eigen_energy * forward_wave_end + side_temp * eigen_energy * the_wave(x_trans_i - 1)) ...
                                             ) ...
                                         ) / (1.0 - side_temp * (the_potential(x_trans_i + 1) - eigen_energy));
                the_wave(1: 1: (x_trans_i - 1)) = the_wave(1: 1: (x_trans_i - 1)) .* (the_wave(x_trans_i + 1) / forward_wave_extension);
                forward_wave_end = forward_wave_end * (the_wave(x_trans_i + 1) / forward_wave_extension);
            end
            was_comparing_deri = true;
            the_diff = ((the_wave(x_trans_i + 1) - the_wave(x_trans_i)) - (forward_wave_end - the_wave(x_trans_i - 1))) / x_step_len ...
                       / sqrt(int_trapezoidal(the_wave .* the_wave, x_step_len));
        else
            was_comparing_deri = false;
            the_diff = (boundary_val - forward_wave_end) / sqrt(int_trapezoidal(the_wave .* the_wave, x_step_len));
        end
        if (not_first_try)
            if (the_diff * first_diff <= 0.0)
                eigen_energy = eigen_energy - e_step_len;
                e_step_len = e_step_len * 0.5;
            end
        else
            first_diff = the_diff;
            not_first_try = true;
        end
        eigen_energy = eigen_energy + e_step_len;
        if (abs(e_step_len) < e_tolerance)
            break;
        end
    end
    is_solution = true;
    if (was_comparing_deri && (abs(the_diff) > deri_diff_tolerance))
        is_solution = false;
        i_state = i_state - 1;
    end
    if (is_solution)
        the_wave = the_wave ./ sqrt(int_simpson(the_wave .* the_wave, x_step_len));
        actual_eigen_val = (round(sqrt(eigen_energy / (pi * pi))) * pi) ^ 2.0;
        if (use_nonzero_potential_eg)
            fprintf("% 1d | %# 13.*f\n", i_state, sig_digits, eigen_energy);
        else
            fprintf("%# 14.*f | %# 13.*f | %# 13.8f%%\n", sig_digits, actual_eigen_val, sig_digits, eigen_energy, 100.0 * (eigen_energy - actual_eigen_val) / actual_eigen_val);
        end
        figure();
        plot(x, the_wave, "k-");
        xlabel("x"); ylabel("u(x)");
        switch (i_state)
            case (1)
                title("基态");
                if (use_nonzero_potential_eg)
                    title("V(x) = 100 * exp(-10 * (x - 0.5)^2), n = 1");
                end
            case (2)
                title("第一激发态");
                if (use_nonzero_potential_eg)
                    title("V(x) = 100 * exp(-10 * (x - 0.5)^2), n = 2");
                end
            case (3)
                title("第二激发态");
                if (use_nonzero_potential_eg)
                    title("V(x) = 100 * exp(-10 * (x - 0.5)^2), n = 3");
                end
            case (4)
                title("第三激发态");
                if (use_nonzero_potential_eg)
                    title("V(x) = 100 * exp(-10 * (x - 0.5)^2), n = 4");
                end
            case (5)
                title("第四激发态");
                if (use_nonzero_potential_eg)
                    title("V(x) = 100 * exp(-10 * (x - 0.5)^2), n = 5");
                end
        end
    end
    if (i_state ~= to_nth_state)
        eigen_energy = eigen_energy + e_step_len;
    end
end

% function integral_val = int_trapezoidal(integrand_func, lower_x, upper_x, step_len)
% function integral_val = int_simpson(integrand_func, lower_x, upper_x, step_len)
