% Three-story shear building modal analysis
% Model: lumped mass 1000 kg per floor, story stiffness 500 kN/m
% SI units throughout

clear; close all;

m = 1000;      % kg per floor
k_story = 500e3; % N/m per story (500 kN/m)

% --- Mass matrix (3 DOFs: floor displacements) ---
M = m * eye(3);

% --- Stiffness matrix for 3-story shear building ---
% K = [k1+k2,  -k2,    0;
%        -k2, k2+k3,  -k3;
%          0,  -k3,   k3]
K = k_story * [2, -1,  0;
              -1,  2, -1;
               0, -1,  1];

% --- Solve generalized eigenvalue problem K*phi = lambda*M*phi ---
[V, D] = eig(K, M);

lambda = diag(D);

% Sort by ascending eigenvalue
[lambda, idx] = sort(lambda);
V = V(:, idx);

% --- Circular frequencies and periods ---
omega = sqrt(lambda);          % rad/s
T = 2 * pi ./ omega;            % s

% --- Normalize mode shapes (max absolute = 1) ---
for i = 1:3
    V(:, i) = V(:, i) / max(abs(V(:, i)));
end

% --- Display results ---
fprintf('=== MATLAB Modal Analysis Results ===\n');
fprintf('Circular frequencies (rad/s):  %.6f  %.6f  %.6f\n', omega);
fprintf('Periods (s):                  %.6f  %.6f  %.6f\n', T);
fprintf('Natural frequencies (Hz):     %.6f  %.6f  %.6f\n', omega/(2*pi));
fprintf('\nNormalized mode shapes (columns = modes):\n');
disp(V);

% --- Write JSON output ---
results.circular_frequencies_rad_per_s = omega';
results.periods_s = T';
results.normalized_mode_shapes = V;

json_str = jsonencode(results);

% Save to the matlab verification directory
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
out_path = fullfile(script_dir, 'modal_results.json');
fid = fopen(out_path, 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

fprintf('\nResults saved to: %s\n', out_path);
