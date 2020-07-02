clear;

% Using a system of reaction-diffusion equations that desbribes a certain pattern-dynamic system.
diffusion_factor = 1.0;
arg_e = 0.03; % In the range 0.01 < arg_e < 0.06, suitable initial conditions lead to steadily rotating spiral waves,
              % and in the range arg_e > 0.07, spiral waves will break up and the system will quickly fall into a turbulence state.
arg_a = 0.84; arg_b = 0.07; % Better do not modify the values of these two arguments.
reaction_x_e_func = @(u, v) (u .* (1.0 - u) .* (u - (v + arg_b) ./ arg_a)); % Reaction function multiplies arg_e.
activator_func = @(u) (double(u > 1.0) + double((u >= 1.0/3.0) & (u <= 1.0)) .* (1.0 - 6.75 .* u .* (u - 1.0).^2.0));

indep_var_count = int32(2); dep_var_count = int32(2);
a = int32(1); b = int32(2); c = int32(3); % Aliases for the boundary constants
boundary_consts = zeros(2 * indep_var_count * dep_var_count, 3); % a * u + b * p/px(u) = c at the boundaries.
u_xl = int32(1); u_xu = int32(2); u_yl = int32(3); u_yu = int32(4); % Aliases: u at lower x, u at upper x, u at lower y, u at upper y.
v_xl = int32(5); v_xu = int32(6); v_yl = int32(7); v_yu = int32(8); % Aliases: v at lower x, v at upper x, v at lower y, v at upper y.
boundary_consts(:, a) = 0.0;
boundary_consts(:, b) = 1.0;
boundary_consts(:, c) = 0.0;

u_init_dist_func = @(x, y) (0.0); % Initial distribution of u.
v_init_dist_func = @(x, y) (0.0); % Initial distribution of v.
nonspecific = int32(0); % Alias: to not use any specific initial distribution below (namely to use the initial distribution above).
spiral_gen = int32(1); % Alias: to use the initial distribution that generates a spiral wave.
the_spiral = int32(2); % Alias: to use the initial distribution that is a spiral wave.
the_turbulence = int32(3); % Alias: to use the initial distribution that is a turbulence state.
use_specific_init = nonspecific;

x_step_len = 0.2;
lower_x = 0.0; upper_x = 50.0;
x = lower_x: x_step_len: upper_x;
x_count = int32(length(x));
if (use_specific_init == spiral_gen)
    x_half = (upper_x + lower_x) * 0.5;
    u_init_dist_func = @(x, y) (1.0 .* double((x > (x_half - 1.0)) & (x < x_half)) .* double(y > x_half));
    v_init_dist_func = @(x, y) (1.0 .* double((x > x_half) & (x < (x_half + 1.0))) .* double(y > x_half));
end

time_step_len = 0.04;
time_length = 40.0;
actual_steps_in_step = int32(8);
actual_t_step_len = time_step_len / double(actual_steps_in_step); % The forward Euler method is stable for diffusion equations when actual_t_step_len < x_step_len ^ 2 / (2 * diffusion_factor).
% In reaction-diffusion equations, effectively, the upper limit of actual_t_step_len might be smaller.
experience_factor_max = 0.25;
d2_term_factor = diffusion_factor * actual_t_step_len / (x_step_len * x_step_len);
if (d2_term_factor >= (experience_factor_max - eps(experience_factor_max)))
    fprintf("WARNING: too large actual_t_step_len (%G), better be smaller than %G.\n", actual_t_step_len, experience_factor_max * x_step_len * x_step_len / diffusion_factor);
end
t_count = int32(time_length / time_step_len) + 1;

u_field_vals = zeros(x_count, x_count, t_count);
v_field_vals = zeros(x_count, x_count, t_count);
for i_y = 1: 1: x_count
    u_field_vals(:, i_y, 1) = u_init_dist_func(x, x(i_y));
    v_field_vals(:, i_y, 1) = v_init_dist_func(x, x(i_y));
end
if (use_specific_init == the_spiral)
    if ((x_step_len == 0.2) && (lower_x == 0.0) && (upper_x == 50.0))
        u_load = load("data_spiral-wave_u(x,y)_e=0.03_x-y-step=0.2_x-y=0~50.mat", "u_last");
        v_load = load("data_spiral-wave_v(x,y)_e=0.03_x-y-step=0.2_x-y=0~50.mat", "v_last");
        % The saved spiral wave data originally evolved with arg_e = 0.03.
        u_field_vals(:, :, 1) = u_load.u_last;
        v_field_vals(:, :, 1) = v_load.v_last;
        clear u_load v_load;
    else
        fprintf("ERROR: the x now does not match the x in the saved spiral wave data.\n");
        return;
    end
elseif (use_specific_init == the_turbulence)
    if ((x_step_len == 0.2) && (lower_x == 0.0) && (upper_x == 50.0))
        u_load = load("data_turbulence_u(x,y)_e=0.12_x-y-step=0.2_x-y=0~50.mat", "u_last");
        v_load = load("data_turbulence_v(x,y)_e=0.12_x-y-step=0.2_x-y=0~50.mat", "v_last");
        % The saved turbulence data originally evolved with arg_e = 0.12.
        u_field_vals(:, :, 1) = u_load.u_last;
        v_field_vals(:, :, 1) = v_load.v_last;
        clear u_load v_load;
    else
        fprintf("ERROR: the x now does not match the x in the saved turbulence data.\n");
        return;
    end
end
u_field_just_now = u_field_vals(:, :, 1);
v_field_just_now = v_field_vals(:, :, 1);

do_progress_printing = true;
prog_section_len = t_count / int32(100);
i_last_section = int32(0);
i_prog_section = int32(0);
if (do_progress_printing)
    fprintf("Calculation starts at time %.2fs.\n", cputime());
end
reaction_factor = actual_t_step_len / arg_e;
for i_time = 2: 1: t_count
    %boundary_consts(u_xl, c) = double(mod(i_time, int32(ceil(2.48 / time_step_len))) <= int32(2.0 / time_step_len));
        % This is a high-frequency travelling wave generator (corroborated with arg_e = 0.03).
        % Do set boundary_consts(u_xl, [a, b]) = [1.0, 0.0] earlier.
    for i_time_inner = 1: 1: actual_steps_in_step
        u_field_vals(2: 1: (end - 1), 2: 1: (end - 1), i_time) = ( ... 
            u_field_just_now(2: 1: (end - 1), 2: 1: (end - 1)) + ( ...
                reaction_factor .* reaction_x_e_func(u_field_just_now(2: 1: (end - 1), 2: 1: (end - 1)), v_field_just_now(2: 1: (end - 1), 2: 1: (end - 1))) ...
                + d2_term_factor .* ( ...
                    (-4.0) .* u_field_just_now(2: 1: (end - 1), 2: 1: (end - 1)) ...
                    + u_field_just_now(2: 1: (end - 1), 3: 1: end) + u_field_just_now(2: 1: (end - 1), 1: 1: (end - 2)) ...
                    + u_field_just_now(3: 1: end, 2: 1: (end - 1)) + u_field_just_now(1: 1: (end - 2), 2: 1: (end - 1)) ...
                ) ...
            ) ...
        );
        v_field_vals(2: 1: (end - 1), 2: 1: (end - 1), i_time) = v_field_just_now(2: 1: (end - 1), 2: 1: (end - 1)) ...
                                                                 + actual_t_step_len .* ( ...
                                                                     activator_func(u_field_just_now(2: 1: (end - 1), 2: 1: (end - 1))) ...
                                                                     - v_field_just_now(2: 1: (end - 1), 2: 1: (end - 1)) ...
                                                                 );
        u_field_vals(1, :, i_time) = (boundary_consts(u_xl, b) .* u_field_vals(2, :, i_time) - boundary_consts(u_xl, c) * x_step_len) ./ (boundary_consts(u_xl, b) - boundary_consts(u_xl, a) * x_step_len);
        u_field_vals(end, :, i_time) = (boundary_consts(u_xu, b) .* u_field_vals(end - 1, :, i_time) + boundary_consts(u_xu, c) * x_step_len) ./ (boundary_consts(u_xu, b) + boundary_consts(u_xu, a) * x_step_len);
        u_field_vals(:, 1, i_time) = (boundary_consts(u_yl, b) .* u_field_vals(:, 2, i_time) - boundary_consts(u_yl, c) * x_step_len) ./ (boundary_consts(u_yl, b) - boundary_consts(u_yl, a) * x_step_len);
        u_field_vals(:, end, i_time) = (boundary_consts(u_yu, b) .* u_field_vals(:, end - 1, i_time) + boundary_consts(u_yu, c) * x_step_len) ./ (boundary_consts(u_yu, b) + boundary_consts(u_yu, a) * x_step_len);
        v_field_vals(1, :, i_time) = (boundary_consts(v_xl, b) .* v_field_vals(2, :, i_time) - boundary_consts(v_xl, c) * x_step_len) ./ (boundary_consts(v_xl, b) - boundary_consts(v_xl, a) * x_step_len);
        v_field_vals(end, :, i_time) = (boundary_consts(v_xu, b) .* v_field_vals(end - 1, :, i_time) + boundary_consts(v_xu, c) * x_step_len) ./ (boundary_consts(v_xu, b) + boundary_consts(v_xu, a) * x_step_len);
        v_field_vals(:, 1, i_time) = (boundary_consts(v_yl, b) .* v_field_vals(:, 2, i_time) - boundary_consts(v_yl, c) * x_step_len) ./ (boundary_consts(v_yl, b) - boundary_consts(v_yl, a) * x_step_len);
        v_field_vals(:, end, i_time) = (boundary_consts(v_yu, b) .* v_field_vals(:, end - 1, i_time) + boundary_consts(v_yu, c) * x_step_len) ./ (boundary_consts(v_yu, b) + boundary_consts(v_yu, a) * x_step_len);
        u_field_just_now = u_field_vals(:, :, i_time);
        v_field_just_now = v_field_vals(:, :, i_time);
    end
    if (do_progress_printing && (i_time ~= t_count))
        i_prog_section = int32(floor(double(i_time) / double(prog_section_len)));
        if (i_prog_section ~= i_last_section)
            if (i_last_section == int32(0))
                fprintf("Calculation in progress. Progress: ##.#%%.\n")
            end
            fprintf("\b\b\b\b\b\b\b%4.1f%%.\n", floor(double(i_time) / double(t_count) * 1000.0) * 0.1);
        end
        i_last_section = i_prog_section;
    end
end
if (do_progress_printing)
    fprintf("\b\b\b\b\b\b\b100.0%%.\nCalculation finishes at time %.2fs.\n", cputime());
end

%%%
lets_record_this = true; % Set this to true to generate avi file.
lets_play_this = false; % Set this to true to watch the diffusion process.
if (lets_record_this || lets_play_this)
    film_scene = figure();
    film_scene.Resize = 'off'; film_scene.WindowStyle = 'modal';
    film_scene.Color = 'white';
    canvas_size = film_scene.InnerPosition;
    canvas_size([3, 4]) = [400.0, 381.0];
    film_scene.InnerPosition = canvas_size;
    the_fps = 1.0 / time_step_len;
    film_mat = moviein(t_count);
    for i_time = 1: 1: t_count
        image(x, x, 256.0 .* u_field_vals(:, :, i_time)');
        axis([lower_x, upper_x, lower_x, upper_x]);
        film_mat(:, i_time) = getframe(film_scene);
    end
    if (lets_play_this)
        movie(film_scene, film_mat, 1, the_fps);
    end
    if (lets_record_this)
        slow_rate = round(2500.0 / the_fps) * 0.01;
        video_cache = VideoWriter("output_film_x" + string(slow_rate) + ".avi");
        video_cache.FrameRate = slow_rate * the_fps;
        open(video_cache);
            writeVideo(video_cache, film_mat);
        close(video_cache);
    end
    close(film_scene);
end
