function eigen_val = symmtridiag_eigen(major_diagonal, side_diagonal, guess_eigen_val, init_step_len, tolerance, max_step)
% Gives an eigen value of a symmetric tridiagonal matrix depending on the choice of the guess eigen value.
% Fu Chunhao, 1 June 2020
%
% function eigen_val = symmtridiag_eigen(major_diagonal, side_diagonal, guess_eigen_val, init_step_len, tolerance, max_step)
%     mat_dim = length(major_diagonal(:)); % Dimension of the tridiagonal matrix.
%     % mat_dim = length(side_diagonal(:)) + 1;
%     step_now = uint32(1); max_step = uint32(max_step);
%     ...
%     while ((max_step == uint32(0)) || (step_now < max_step))
%         step_now = step_now + 1;
%         ... % (break by meeting the tolerance)
    mat_dim = length(major_diagonal(:)); % Dimension of the tridiagonal matrix.
    if (int32(mat_dim) == int32(1))
        eigen_val = major_diagonal(1);
        return;
    end
    step_now = uint32(1); max_step = uint32(max_step);
    eigen_val = guess_eigen_val;
    step_len = init_step_len;
    poly_1 = major_diagonal(1) - eigen_val;
    poly_2 = (major_diagonal(2) - eigen_val) * poly_1 - side_diagonal(1) * side_diagonal(1);
    for i_order = 3: 1: mat_dim
        poly_3 = (major_diagonal(i_order) - eigen_val) * poly_2 - side_diagonal(i_order - 1) * side_diagonal(i_order - 1) * poly_1;
        poly_1 = poly_2;
        poly_2 = poly_3;
    end
    init_poly_val = poly_2;
    while ((max_step == uint32(0)) || (step_now < max_step))
        step_now = step_now + 1;
        if (abs(step_len) < tolerance)
            break;
        end
        eigen_val = eigen_val + step_len;
        poly_1 = major_diagonal(1) - eigen_val;
        poly_2 = (major_diagonal(2) - eigen_val) * poly_1 - side_diagonal(1) * side_diagonal(1);
        for i_order = 3: 1: mat_dim
            poly_3 = (major_diagonal(i_order) - eigen_val) * poly_2 - side_diagonal(i_order - 1) * side_diagonal(i_order - 1) * poly_1;
            poly_1 = poly_2;
            poly_2 = poly_3;
        end
        if (init_poly_val * poly_2 < 0.0)
            eigen_val = eigen_val - step_len;
            step_len = 0.5 * step_len;
        end
    end
end
