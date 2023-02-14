clear
clc
close all

%% Define specifications
% Calculate number of stages
inlet_temp_F = 90;
outlet_temp_F = 105;
num_of_stages = 6;

% Temperature (tear variable)
% Assumed a linear variation of temperature with each stage
temp_K = linspace(inlet_temp_F, outlet_temp_F, num_of_stages)';

% Vapour flowrate (tear variable)
vapour_rate_lbmolh = zeros(num_of_stages, 1);
vapour_rate_lbmolh(:) = 800; % lbmol/h

% Stage pressure
% Assumed no pressure drop
pressure_psia = zeros(num_of_stages, 1);
pressure_psia(:) = 400;

% Feed molar compisition
zC1 = 160/800;
zC2 = 370/800;
zC3 = 240/800;
zC4 = 25/800;
zC5 = 5/800;
zOil = 0;

feed_mol_comp = [zC1 zC2 zC3 zC4 zC5 zOil];

comp_number = length(feed_mol_comp);

% Feed molar flow rate
feed_rate_lbmolh = zeros(num_of_stages, 1);
feed_rate_lbmolh(1) = 165;
feed_rate_lbmolh(6) = 800; % lbmol/h

% Create matrix for feed composition for each stage (useful for enthalpy
% calculation)
feed_mol_comp_matrix = zeros(num_of_stages, comp_number);
feed_mol_comp_matrix(1, :) = [0 0 0 0.05/165 0.78/165 164.17/165];
feed_mol_comp_matrix(6, :) = feed_mol_comp;

% Liquid side stream
U = zeros(num_of_stages, 1);

% Equilibrium constant
eq_const = zeros(num_of_stages, comp_number);
for j = 1:num_of_stages

    eq_const(j, :) = [6.65 1.64 0.584 0.195 0.0713 0.0001];

end

%% Burningham-Otto algorithm loop
max_iteration = 10;

for iter = 1:max_iteration
%%%%%%%%%%%%%%%%%%%%%% Mass balance equation solving %%%%%%%%%%%%%%%%%%%%%%
% Calculate x matrix coefficients
[A_x_matrix, B_x_matrix, C_x_matrix, D_x_matrix] = xMatrixCoefficientsExample(num_of_stages, comp_number, vapour_rate_lbmolh, feed_rate_lbmolh,...
                                                    U, eq_const, feed_mol_comp_matrix);

% Solve for x values using Thomas function
liquid_comp = zeros(num_of_stages, comp_number);

for i = 1:comp_number

    liquid_comp(:, i) = tridiagonal_vector(A_x_matrix, B_x_matrix(:, i), C_x_matrix(:, i), D_x_matrix(:, i));

end

% Liquid flow rate calculated using vapour flow rate
liquid_rate_lbmolh = zeros(num_of_stages, 1);

for j = 1:num_of_stages - 1 % -1 to prevent over-indexing due to V(j + 1)

    liquid_rate_lbmolh(j) = vapour_rate_lbmolh(j + 1) + (sum(feed_rate_lbmolh(1:j)) - sum(U(1:j)))...
                            - vapour_rate_lbmolh(1);

end

% Calculate last element of liquid flow rate matrix
j = num_of_stages;
liquid_rate_lbmolh(j) = (sum(feed_rate_lbmolh(1:j)) - sum(U(1:j))) - vapour_rate_lbmolh(1);

% Calculate adjusted liquid flow rate using sum-rates equation
for j = 1:num_of_stages

    liquid_rate_lbmolh(j) = liquid_rate_lbmolh(j)*sum(liquid_comp(j, :));

end

% Calculate adjusted vapour flow rate
for j = 2:num_of_stages

    vapour_rate_lbmolh(j) = liquid_rate_lbmolh(j - 1) - liquid_rate_lbmolh(num_of_stages) + (sum(feed_rate_lbmolh(j:num_of_stages))...
                            - sum(U(j:num_of_stages)));

end

% Calculate first element of vapour flow rate (L(j - 1) = 0)
j = 1;
vapour_rate_lbmolh(j) = - liquid_rate_lbmolh(num_of_stages) + (sum(feed_rate_lbmolh(j:num_of_stages))...
                        - sum(U(j:num_of_stages)));

% Normalise x values
for i = 1:comp_number
    for j = 1:num_of_stages

        liquid_comp(j, i) = liquid_comp(j, i)/sum(liquid_comp(j, :));

    end
end

% Calculate y values
vapour_comp = zeros(num_of_stages, comp_number);

for i = 1:comp_number
    for j = 1:num_of_stages

        vapour_comp(j, i) = eq_const(j, i)*liquid_comp(j, i);

    end
end

% %%%%%%%%%%%%%%%% double check if you normalise y first THEN calculate K

% Normalise y values
for i = 1:comp_number
    for j = 1:num_of_stages

        vapour_comp(j, i) = vapour_comp(j, i)/sum(vapour_comp(j, :));

    end
end

% % Calculate adjusted equilibrium constant values
% for i = 1:comp_number
%     for j = 1:num_of_stages
% 
%         eq_const(j, i) = vapour_comp(j, i)/liquid_comp(j, i);
%         
%     end
% end
% 
end
% %%%%%%%%%%%%%%%%%%%%%% Heat balance equation solving %%%%%%%%%%%%%%%%%%%%%%
% % Define Cp coefficients for liquid and vapour components (Cp = A + BT ...)
% % Values obtained from NIST (and polyfit data), DDBST, Kasprzycka-Guttman &
% % Odzeniak (1991)
% % All in J/mol K
% cp_coeff_water_liquid = [56.9866260269178 0.187904654544170 -0.000654742631656604 7.70428331682889e-07]';
% cp_coeff_water_vapour = [33.9369221307440 -0.00797039884420184 2.64722202554198e-05 -1.12847969807767e-08]';
% 
% cp_coeff_hexane_liquid = [172.119999999999 -0.183779999999992 0.000887339999999974 2.46146849386985e-20]';
% cp_coeff_hexane_vapour = [-4.44257770284955 0.587668036942747 -0.000321385755768740 6.86338615245967e-08]';
% 
% cp_coeff_oil = [14528.8966697806 -102.025824787746 0.268637783393027 -0.000230599204716798]';
% 
% cp_coeff_meal = [meal_cp_coefficient(mmMeal)]'; % Used a function because unsure of what mmMeal is
% 
% % Create matrix of Cp coefficients
% cp_coeff_liquid = [cp_coeff_oil cp_coeff_water_liquid cp_coeff_meal cp_coeff_hexane_liquid];
% cp_coeff_vapour = zeros(length(cp_coeff_water_liquid), comp_number);
% cp_coeff_vapour(:, 2) = cp_coeff_water_vapour;
% cp_coeff_vapour(:, 4) = cp_coeff_hexane_vapour;
% 
% % Calculate vapour and liquid enthalpies
% vapour_enthalpy = zeros(num_of_stages, 1);
% liquid_enthalpy = zeros(num_of_stages, 1);
% feed_enthalpy = zeros(num_of_stages, 1);
% directSteam_enthalpy = zeros(num_of_stages, 1);
% coeff_sum_vapour = 0;
% coeff_sum_liquid = 0;
% RHS_vapour = 0;
% RHS_liquid = 0;
% RHS_feed = 0;
% RHS_directSteam = 0;
% 
% for j = 1:num_of_stages
%     for i = 1:comp_number
%         for m = 1:length(cp_coeff_water_liquid)
%             
%             % Calculates the A(i) + B(i)T + C(i)T^2 ... function
%             coeff_sum_vapour = coeff_sum_vapour + (cp_coeff_vapour(m, i)*temp_K(j)^(m - 1)); 
%             coeff_sum_liquid = coeff_sum_liquid + (cp_coeff_liquid(m, i)*temp_K(j)^(m - 1));
%                 
%         end
%         
%         % Calculates the RHS function for 1 component and sums it to the previous component
%         RHS_vapour = RHS_vapour + vapour_comp(j, i)*coeff_sum_vapour; 
%         RHS_liquid = RHS_liquid + liquid_comp(j, i)*coeff_sum_liquid;
%         RHS_feed = RHS_feed + feed_mol_comp_matrix(j, i)*coeff_sum_liquid;
%         RHS_directSteam = RHS_directSteam + directSteam_comp(j, i)*coeff_sum_vapour;
%         
%         % Resets to 0 for use in next component
%         coeff_sum_vapour = 0;
%         coeff_sum_liquid = 0;
% 
%     end
% 
%     vapour_enthalpy(j) = RHS_vapour;
%     liquid_enthalpy(j) = RHS_liquid;
%     feed_enthalpy(j) = RHS_feed;
%     directSteam_enthalpy(j) = RHS_directSteam;
% 
%     RHS_vapour = 0;
%     RHS_liquid = 0;
%     RHS_feed = 0;
%     RHS_directSteam = 0;
% 
% end
% 
% % Calculate vapour and liquid enthalpy partial derivatives (dH/dT)
% vapour_enthalpy_derivative = zeros(num_of_stages, 1);
% liquid_enthalpy_derivative = zeros(num_of_stages, 1);
% coeff_sum_vapour = 0;
% coeff_sum_liquid = 0;
% RHS_vapour = 0;
% RHS_liquid = 0;
% 
% for j = 1:num_of_stages
%     for i = 1:comp_number
%         for m = 2:length(cp_coeff_water_liquid)
%             
%             % Calculates the B(i) + 2C(i)T ... function
%             coeff_sum_vapour = coeff_sum_vapour + ((m - 1)*cp_coeff_vapour(m, i)*temp_K(j)^(m - 2)); 
%             coeff_sum_liquid = coeff_sum_liquid + ((m - 1)*cp_coeff_liquid(m, i)*temp_K(j)^(m - 2));
%                 
%         end
%         
%         % Calculates the RHS function for 1 component and sums it to the previous component
%         RHS_vapour = RHS_vapour + vapour_comp(j, i)*coeff_sum_vapour; 
%         RHS_liquid = RHS_liquid + liquid_comp(j, i)*coeff_sum_liquid; 
%         
%         % Resets to 0 for use in next component
%         coeff_sum_vapour = 0;
%         coeff_sum_liquid = 0;
% 
%     end
% 
%     vapour_enthalpy_derivative(j) = RHS_vapour;
%     liquid_enthalpy_derivative(j) = RHS_liquid;
% 
%     RHS_vapour = 0;
%     RHS_liquid = 0;
% 
% end
% 
% % Calculate temperature tridiagonal matrix coefficients
% A_T_vector = zeros(num_of_stages - 1, 1);
% B_T_vector = zeros(num_of_stages, 1);
% C_T_vector = zeros(num_of_stages - 1, 1);
% D_T_vector = zeros(num_of_stages, 1);
% 
% % Populate A vector
% for j = 2:num_of_stages
% 
%     A_T_vector(j - 1) = liquid_rate_molsec(j)*liquid_enthalpy_derivative(j);
% 
% end
% 
% % Populate B vector
% for j = 1:num_of_stages
% 
%     B_T_vector(j) = -liquid_rate_molsec(j)*liquid_enthalpy_derivative(j) ...
%                     -(vapour_rate_molsec(j) + directSteam_rate_molsec(j))*vapour_enthalpy_derivative(j);
% 
% end
% 
% % Populate C vector
% for j = 1:num_of_stages - 1
% 
%     C_T_vector(j) = vapour_rate_molsec(j + 1)*vapour_enthalpy_derivative(j + 1);
% 
% end
% 
% % Populate D vector
% for j = 2:num_of_stages - 1
% 
%     D_T_vector(j) = liquid_rate_molsec(j - 1)*liquid_enthalpy(j - 1) + vapour_rate_molsec(j + 1)*vapour_enthalpy(j + 1)...
%                     + feed_rate_molsec(j)*feed_enthalpy(j) + directSteam_rate_molsec(j)*directSteam_enthalpy(j)...
%                     + indirectSteam_energy(j) - liquid_rate_molsec(j)*liquid_enthalpy(j) - vapour_rate_molsec(j)*vapour_enthalpy(j);
% 
% end
% 
% % Calculate D(1), where L(j - 1) = 0;
% j = 1;
% D_T_vector(1) = vapour_rate_molsec(j + 1)*vapour_enthalpy(j + 1)...
%                 + feed_rate_molsec(j)*feed_enthalpy(j) + directSteam_rate_molsec(j)*directSteam_enthalpy(j)...
%                 + indirectSteam_energy(j) - liquid_rate_molsec(j)*liquid_enthalpy(j) - vapour_rate_molsec(j)*vapour_enthalpy(j);
% 
% % Calculate D(end), where V(j + 1) = 0;
% j = length(D_T_vector);
% D_T_vector(j) = liquid_rate_molsec(j - 1)*liquid_enthalpy(j - 1) + feed_rate_molsec(j)*feed_enthalpy(j)...
%                 + directSteam_rate_molsec(j)*directSteam_enthalpy(j) + indirectSteam_energy(j)...
%                 - liquid_rate_molsec(j)*liquid_enthalpy(j) - vapour_rate_molsec(j)*vapour_enthalpy(j);
% 
% % Solve temperature tridiagonal matrix using Thomas function
% deltaT = tridiagonal_vector(A_T_vector, B_T_vector, C_T_vector, D_T_vector);
% 
% % Add temperature correction to stage temperature
% temp_K = temp_K + deltaT;
% 
% % Calculate exit criteria
% exit_criteria_tau = sum(deltaT.^2);
% exit_criteria_N = 0.01*num_of_stages;
% 
% % Check if exit criteria is met
% if exit_criteria_tau <= exit_criteria_N
% 
%     break
% 
% end
% 
% % Check temperature matrix coefficients
% 
% end