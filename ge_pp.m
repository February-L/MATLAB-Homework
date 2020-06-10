function x = ge_pp(augmented_matrix, x_count)
% Solves a system of linear equations by using the Gaussian elimination with partial pivoting.
% Fu Chunhao, 23 April 2020
%
% function x = ge_pp(augmented_matrix, x_count)
%     aug_mat = augmented_matrix(1: 1: x_count, 1: 1: (x_count + 1));
%     x = zeros(x_count, 1);
    aug_mat = augmented_matrix(1: 1: x_count, 1: 1: (x_count + 1));
    x = zeros(x_count, 1);
    for i_col = 1: 1: (x_count - 1)
        i_row_prcpl = i_col;
        for i_row = (i_col + 1): 1: x_count
            if (abs(aug_mat(i_row, i_col)) > abs(aug_mat(i_row_prcpl, i_col)))
                i_row_prcpl = i_row;
            end
        end
        row_temp = aug_mat(i_col, :);
        aug_mat(i_col, :) = aug_mat(i_row_prcpl, :);
        aug_mat(i_row_prcpl, :) = row_temp;
        aug_mat((i_col + 1): 1: x_count, :) = aug_mat((i_col + 1): 1: x_count, :) - aug_mat((i_col + 1): 1: x_count, i_col) * aug_mat(i_col, :) ./ aug_mat(i_col, i_col) ;
    end
    x(x_count, 1) = aug_mat(x_count, x_count + 1) / aug_mat(x_count, x_count);
    for i_row = (x_count - 1): (-1): 1
        x(i_row, 1) = (aug_mat(i_row, x_count + 1) - aug_mat(i_row, (i_row + 1): 1: x_count) * x((i_row + 1): 1: x_count, 1)) / aug_mat(i_row, i_row);
    end
end
