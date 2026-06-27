%% Modal analysis of a 3-story shear building (MATLAB)
% Each story: lumped mass m = 1000 kg, story stiffness k = 500 kN/m
% SI units throughout.

clear; clc;

m = 1000;          % kg per floor
k = 500e3;         % N/m per story

%% Assemble mass and stiffness matrices
M = m * eye(3);

K = k * [ 2, -1,  0;
         -1,  2, -1;
          0, -1,  1];

%% Solve generalized eigenvalue problem  K*phi = lambda*M*phi
[Phi, Lambda] = eig(K, M);

omega2 = diag(Lambda);                 % eigenvalues (rad/s)^2
omega  = sqrt(omega2);                 % circular frequencies rad/s
T      = 2*pi ./ omega;                % periods s

% Sort by ascending frequency
[omega, idx] = sort(omega);
T = T(idx);
Phi = Phi(:, idx);

% Normalize mode shapes so that top-floor amplitude = 1
for j = 1:3
    Phi(:, j) = Phi(:, j) / Phi(3, j);
end

%% Display results
fprintf('=== MATLAB Modal Analysis Results ===\n\n');
for j = 1:3
    fprintf('Mode %d:\n', j);
    fprintf('  omega  = %.4f rad/s\n', omega(j));
    fprintf('  f      = %.4f Hz\n', omega(j)/(2*pi));
    fprintf('  T      = %.6f s\n', T(j));
    fprintf('  phi    = [%.6f, %.6f, %.6f]\n\n', Phi(1,j), Phi(2,j), Phi(3,j));
end

%% Export to JSON
json_struct.circular_frequencies_rad_per_s = omega';
json_struct.periods_s = T';
json_struct.normalized_mode_shapes = Phi';

json_str = jsonencode(json_struct);

fid = fopen('modal_results.json', 'w');
if fid == -1
    error('Cannot open modal_results.json for writing');
end
fprintf(fid, '%s', json_str);
fclose(fid);

fprintf('Results written to modal_results.json\n');
