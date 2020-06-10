function x = lu_fact(augmented_matrix, x_count)
% Solves a system of linear equations by using the LU factorisation.
% Fu Chunhao, 23 April 2020
%
% function x = lu_fact(augmented_matrix, x_count)
%     mat_L = zeros(x_count, x_count); % Lower triangular matrix L.
%     mat_U = zeros(x_count, x_count); % Upper triangular matrix U.
%     y = zeros(x_count, 1);
%     x = zeros(x_count, 1);
    mat_L = zeros(x_count, x_count); % Lower triangular matrix L.
    mat_U = zeros(x_count, x_count); % Upper triangular matrix U.
    y = zeros(x_count, 1);
    x = zeros(x_count, 1);
    mat_U(1, :) = augmented_matrix(1, 1: 1: x_count);
    mat_L(2: 1: x_count, 1) = augmented_matrix(2: 1: x_count, 1) ./ mat_U(1, 1);
    mat_L(1, 1) = 1.0;
    for i_elem = 2: 1: x_count
        mat_L(i_elem, i_elem) = 1.0;
        mat_U(i_elem, i_elem: 1: x_count) = augmented_matrix(i_elem, i_elem: 1: x_count) - ...
                                            mat_L(i_elem, 1: 1: (i_elem - 1)) * mat_U(1: 1: (i_elem - 1), i_elem: 1: x_count);
        mat_L((i_elem + 1): 1: x_count, i_elem) = (augmented_matrix((i_elem + 1): 1: x_count, i_elem) - ...
                                                  mat_L((i_elem + 1): 1: x_count, 1: 1: (i_elem - 1)) * mat_U(1: 1: (i_elem - 1), i_elem)) ...
                                                  ./ mat_U(i_elem, i_elem);
    end
    y(1, 1) = augmented_matrix(1, x_count + 1);
    for i_elem = 2: 1: x_count
        y(i_elem, 1) = augmented_matrix(i_elem, x_count + 1) - mat_L(i_elem, 1: 1: (i_elem - 1)) * y(1: 1: (i_elem - 1), 1);
    end
    x(x_count, 1) = y(x_count, 1) / mat_U(x_count, x_count);
    for i_elem = (x_count - 1): (-1): 1
        x(i_elem, 1) = (y(i_elem, 1) - mat_U(i_elem, (i_elem + 1): 1: x_count) * x((i_elem + 1): 1: x_count, 1)) ./ mat_U(i_elem, i_elem);
    end
end
