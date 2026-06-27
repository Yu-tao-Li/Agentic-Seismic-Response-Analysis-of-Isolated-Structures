% Modal analysis of a 3-story shear building
% Each story: lumped mass m = 1000 kg, story stiffness k = 500 kN/m = 500000 N/m

clear; clc;

m = 1000;       % kg per floor
k = 500e3;      % N/m per story

% Mass matrix (lumped)
M = m * eye(3);

% Stiffness matrix for shear building
K = k * [ 2, -1,  0;
         -1,  2, -1;
          0, -1,  1];

% Generalized eigenvalue problem: K*phi = lambda*M*phi
[phi, Lambda] = eig(K, M);

% Extract eigenvalues
omega2 = diag(Lambda);
omega = sqrt(omega2);           % circular frequencies (rad/s)
T = 2 * pi ./ omega;            % periods (s)

% Sort by ascending frequency
[omega, idx] = sort(omega);
T = T(idx);
phi = phi(:, idx);

% Normalize mode shapes (each column to unit max absolute value)
for j = 1:3
    phi(:, j) = phi(:, j) / max(abs(phi(:, j)));
end

% Display results
fprintf('=== 3-Story Shear Building Modal Analysis (MATLAB) ===\n\n');
for j = 1:3
    fprintf('Mode %d:\n', j);
    fprintf('  omega = %.6f rad/s\n', omega(j));
    fprintf('  T     = %.6f s\n', T(j));
    fprintf('  phi   = [%.6f, %.6f, %.6f]\n\n', phi(1,j), phi(2,j), phi(3,j));
end

% Write JSON results
json.frequencies = omega';
json.periods = T';
json.mode_shapes = phi;

json_str = '{\n';
json_str = [json_str, sprintf('  "circular_frequencies_rad_per_s": [%.10f, %.10f, %.10f],\n', omega(1), omega(2), omega(3))];
json_str = [json_str, sprintf('  "periods_s": [%.10f, %.10f, %.10f],\n', T(1), T(2), T(3))];
json_str = [json_str, '  "normalized_mode_shapes": [\n'];
for j = 1:3
    comma = '';
    if j < 3, comma = ','; end
    json_str = [json_str, sprintf('    [%.10f, %.10f, %.10f]%s\n', phi(1,j), phi(2,j), phi(3,j), comma)];
end
json_str = [json_str, '  ]\n}'];

fid = fopen('modal_results.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

fprintf('Results saved to modal_results.json\n');
