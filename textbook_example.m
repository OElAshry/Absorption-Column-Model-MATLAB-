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

% Normalise y values
for i = 1:comp_number
    for j = 1:num_of_stages

        vapour_comp(j, i) = vapour_comp(j, i)/sum(vapour_comp(j, :));

    end
end

end
