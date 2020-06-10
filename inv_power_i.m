function [eigenvector, eigenvalue] = inv_power_i(the_matrix, guess_eigenvector, guess_eigenvalue, eigenvalue_tolerance, max_step)
% Gives the given matrix's nearest eigen-value to the guess and its corresponding (normalised) eigen-vector by using the inverse power iteration.
% Note that the inverse power iteration method requires (the_matrix - guess_eigenvalue * eye(size(the_matrix))) to be invertible.
% Fu Chunhao, 1 June 2020
%
% function [eigenvector, eigenvalue] = inv_power_i(the_matrix, guess_eigenvector, guess_eigenvalue, eigenvalue_tolerance, max_step)
%     % guess_eigenvector should be a column vector.
%     step_now = uint32(1); max_step = uint32(max_step);
%     the_matrix = the_matrix - guess_eigenvalue * eye(size(the_matrix));
%     ...
%     while ((max_step == uint32(0)) || (step_now < max_step))
%         step_now = step_now + 1;
%         ... % (break by meeting the tolerance)
    step_now = uint32(1); max_step = uint32(max_step);
    the_matrix = the_matrix - guess_eigenvalue * eye(size(the_matrix));
    eigenvector = guess_eigenvector ./ norm(guess_eigenvector);
    vector_temp = the_matrix \ eigenvector;
    last_eigenval = 1.0 / (eigenvector' * vector_temp) + guess_eigenvalue;
    eigenvalue = last_eigenval;
    while ((max_step == uint32(0)) || (step_now < max_step))
        step_now = step_now + 1;
        eigenvector = vector_temp ./ norm(vector_temp);
        vector_temp = the_matrix \ eigenvector;
        eigenvalue = 1.0 / (eigenvector' * vector_temp) + guess_eigenvalue;
        if (abs(eigenvalue - last_eigenval) < eigenvalue_tolerance)
            break;
        end
        last_eigenval = eigenvalue;
    end
    eigenvector = vector_temp ./ norm(vector_temp);
end
