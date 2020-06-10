function x = tridiag_chasing(lower_diagonal, major_diagonal, upper_diagonal, inhomogeneous)
% Solves a system of linear equations whose coeffecient matrix is tridiagonal by using the chasing method.
% Fu Chunhao, 25 April 2020
%
% function x = tridiag_chasing(lower_diagonal, major_diagonal, upper_diagonal, inhomogeneous)
%     x_count = length(major_diagonal(:));
%     % x_count = length(lower_diagonal(:)) + 1; x_count = length(upper_diagonal(:)) + 1;
%     % x_count = length(inhomogeneous(:));
%     x = zeros(x_count, 1);
    x_count = length(major_diagonal(:));
    x = zeros(x_count, 1);
    if (int32(x_count) == int32(1))
        x(1, 1) = inhomogeneous(1) / major_diagonal(1);
        return;
    end
    for i = 1: 1: (x_count - 1)
        major_diagonal(i + 1) = major_diagonal(i + 1) - lower_diagonal(i) * upper_diagonal(i) / major_diagonal(i);
        inhomogeneous(i + 1) = inhomogeneous(i + 1) - lower_diagonal(i) * inhomogeneous(i) / major_diagonal(i);
    end
    x(x_count, 1) = inhomogeneous(x_count) / major_diagonal(x_count);
    for i = (x_count - 1): (-1): 1
        x(i, 1) = (inhomogeneous(i) - upper_diagonal(i) * x(i + 1, 1)) / major_diagonal(i);
    end
end
