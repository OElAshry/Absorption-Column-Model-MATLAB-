function [A_x_matrix, B_x_matrix, C_x_matrix, D_x_matrix] = xMatrixCoefficientsExample(num_of_stages, comp_number, vapour_rate_lbmolh, feed_rate_lbmolh,...
    U, eq_const, feed_mol_comp_matrix)

A_x_matrix = zeros(num_of_stages - 1, 1);
B_x_matrix = zeros(num_of_stages, comp_number);
C_x_matrix = zeros(num_of_stages - 1, comp_number);
D_x_matrix = zeros(num_of_stages, comp_number);

% Populate A matrix
for j = 2:num_of_stages

    A_x_matrix(j - 1) = vapour_rate_lbmolh(j) + (sum(feed_rate_lbmolh(1:j - 1)) - sum(U(1:j - 1)))...
        - vapour_rate_lbmolh(1);

end
   
% Populate B matrix
for i = 1:comp_number
    for j = 1:num_of_stages - 1 % -1 to prevent over-indexing due to V(j + 1)
    
        B_x_matrix(j, i) = -(vapour_rate_lbmolh(j + 1) + (sum(feed_rate_lbmolh(1:j)) - sum(U(1:j)))...
            - vapour_rate_lbmolh(1) + U(j) + eq_const(j, i)*(vapour_rate_lbmolh(j)));
                                                                                           
                                                                                           % Added the above term (based on general function in textbook)
    
    end

    % Calculate last element of the B matrix
    j = num_of_stages;
    B_x_matrix(j, i) = -((sum(feed_rate_lbmolh(1:j)) - sum(U(1:j)))...
        - vapour_rate_lbmolh(1) + U(j) + eq_const(j, i)*(vapour_rate_lbmolh(j)));
end

% Populate C matrix
for i = 1:comp_number
    for j = 1:num_of_stages - 1

        C_x_matrix(j, i) = eq_const(j + 1, i)*vapour_rate_lbmolh(j + 1);

    end
end

% Populate D matrix
for i = 1:comp_number
    for j = 1:num_of_stages

        D_x_matrix(j, i) = -feed_rate_lbmolh(j)*feed_mol_comp_matrix(j, i);
    end
end

end