clear; clc;

script_dir = fileparts(mfilename('fullpath'));
root_dir = fullfile(script_dir, '..', '..');
gm_path = fullfile(root_dir, 'input_data', 'GM1.txt');

dt = 0.02;
g = 9.81;
target_pga = 0.40 * g;
mass_per_story = 1000.0;
story_stiffness = 500.0e3;
zeta = 0.05;

raw_gm = load(gm_path);
raw_gm = raw_gm(:);
n = length(raw_gm);
time = (0:n-1)' * dt;
raw_peak = max(abs(raw_gm));
if raw_peak <= 0
    error('GM1.txt has zero peak acceleration.');
end
ag = raw_gm / raw_peak * target_pga;
input_pga = max(abs(ag));

M = mass_per_story * eye(3);
k = story_stiffness;
K = [ 2*k, -k,    0;
     -k,  2*k, -k;
      0,  -k,   k];

[~, D] = eig(K, M);
omega = sort(sqrt(diag(D)));
rayleigh_matrix = [1/(2*omega(1)), omega(1)/2;
                   1/(2*omega(2)), omega(2)/2];
rayleigh = rayleigh_matrix \ [zeta; zeta];
alpha_m = rayleigh(1);
beta_k = rayleigh(2);
C = alpha_m * M + beta_k * K;

gamma = 0.5;
beta = 0.25;
u = zeros(3, n);
v = zeros(3, n);
a_rel = zeros(3, n);
r = ones(3, 1);
p = -M * r * ag';

a_rel(:, 1) = M \ (p(:, 1) - C * v(:, 1) - K * u(:, 1));

a0 = 1.0 / (beta * dt^2);
a1 = gamma / (beta * dt);
a2 = 1.0 / (beta * dt);
a3 = 1.0 / (2.0 * beta) - 1.0;
a4 = gamma / beta - 1.0;
a5 = dt * (gamma / (2.0 * beta) - 1.0);
Keff = K + a0 * M + a1 * C;

for i = 1:n-1
    peff = p(:, i+1) ...
        + M * (a0 * u(:, i) + a2 * v(:, i) + a3 * a_rel(:, i)) ...
        + C * (a1 * u(:, i) + a4 * v(:, i) + a5 * a_rel(:, i));
    u(:, i+1) = Keff \ peff;
    a_rel(:, i+1) = a0 * (u(:, i+1) - u(:, i)) - a2 * v(:, i) - a3 * a_rel(:, i);
    v(:, i+1) = v(:, i) + dt * ((1.0 - gamma) * a_rel(:, i) + gamma * a_rel(:, i+1));
end

top_disp = u(3, :)';
top_acc_abs = a_rel(3, :)' + ag;

writematrix([time, top_disp], fullfile(script_dir, 'topStoDis.txt'), 'Delimiter', 'tab');
writematrix([time, top_acc_abs], fullfile(script_dir, 'topStoAcc.txt'), 'Delimiter', 'tab');
writematrix([time, ag], fullfile(script_dir, 'normalized_ground_acceleration.txt'), 'Delimiter', 'tab');

response = struct();
response.software = 'MATLAB';
response.model = 'three_story_linear_shear_building_relative_coordinates';
response.output_definition = struct( ...
    'top_story_displacement', 'relative displacement of story 3 with respect to ground, m', ...
    'top_story_acceleration', 'absolute acceleration of story 3, m/s^2');
response.time = time';
response.top_story_displacement = top_disp';
response.top_story_acceleration = top_acc_abs';
response.normalized_ground_acceleration = ag';
response.input_pga_after_normalization = input_pga;
response.dt = dt;
response.num_samples = n;
response.mass_per_story_kg = mass_per_story;
response.story_stiffness_N_per_m = story_stiffness;
response.damping_ratio = zeta;
response.rayleigh_alpha_m = alpha_m;
response.rayleigh_beta_k = beta_k;
response.circular_frequencies_rad_s = omega';
response.newmark_gamma = gamma;
response.newmark_beta = beta;

fid = fopen(fullfile(script_dir, 'response.json'), 'w');
fprintf(fid, '%s', jsonencode(response));
fclose(fid);

fprintf('MATLAB analysis complete. Samples: %d, dt: %.8f, normalized PGA: %.8f m/s^2\n', n, dt, input_pga);
fprintf('Rayleigh alpha_m: %.12g, beta_k: %.12g\n', alpha_m, beta_k);
fprintf('Peak |top displacement|: %.12g m\n', max(abs(top_disp)));
fprintf('Peak |top absolute acceleration|: %.12g m/s^2\n', max(abs(top_acc_abs)));
