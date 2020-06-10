function [eigenvector, eigenvalue] = power_i(the_matrix, guess_eigenvector, eigenvalue_tolerance, max_step)
% Gives the dominant eigen-value and (normalised) eigen-vector of the given matrix by using the power iteration.
% Fu Chunhao, 1 June 2020
%
% function [eigenvector, eigenvalue] = power_i(the_matrix, guess_eigenvector, eigenvalue_tolerance, max_step)
%     % guess_eigenvector should be a column vector.
%     step_now = uint32(1); max_step = uint32(max_step);
%     ...
%     while ((max_step == uint32(0)) || (step_now < max_step))
%         step_now = step_now + 1;
%         ... % (break by meeting the tolerance)
    step_now = uint32(1); max_step = uint32(max_step);
    eigenvector = guess_eigenvector ./ norm(guess_eigenvector);
    vector_temp = the_matrix * eigenvector;
    last_eigenval = eigenvector' * vector_temp;
    eigenvalue = last_eigenval;
    while ((max_step == uint32(0)) || (step_now < max_step))
        step_now = step_now + 1;
        eigenvector = vector_temp ./ norm(vector_temp);
        vector_temp = the_matrix * eigenvector;
        eigenvalue = eigenvector' * vector_temp;
        if (abs(eigenvalue - last_eigenval) < eigenvalue_tolerance)
            break;
        end
        last_eigenval = eigenvalue;
    end
    eigenvector = vector_temp ./ norm(vector_temp);
end
