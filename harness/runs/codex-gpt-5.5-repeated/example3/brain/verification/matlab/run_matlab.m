clear;
clc;

script_dir = fileparts(mfilename('fullpath'));
input_path = fullfile(script_dir, '..', '..', 'input_data', 'GM1.txt');

dt = 0.02;
g = 9.80665;
target_pga = 0.40 * g;
zeta = 0.05;

raw_gm = readmatrix(input_path);
raw_gm = raw_gm(:);
raw_peak = max(abs(raw_gm));
if raw_peak <= 0
    error('GM1.txt has zero peak acceleration.');
end

ag = raw_gm / raw_peak * target_pga;
n = numel(ag);
time = (0:n-1)' * dt;

masses = [500; 1000; 1000; 1000];
stiffnesses = [200e3; 500e3; 500e3; 500e3];
c_iso = 2e3;

M = diag(masses);
K = [stiffnesses(1) + stiffnesses(2), -stiffnesses(2), 0, 0; ...
    -stiffnesses(2), stiffnesses(2) + stiffnesses(3), -stiffnesses(3), 0; ...
    0, -stiffnesses(3), stiffnesses(3) + stiffnesses(4), -stiffnesses(4); ...
    0, 0, -stiffnesses(4), stiffnesses(4)];

[~, D] = eig(K, M);
omega = sort(sqrt(diag(D)));
rayleigh_matrix = [1/(2*omega(1)), omega(1)/2; 1/(2*omega(2)), omega(2)/2];
rayleigh = rayleigh_matrix \ [zeta; zeta];
alpha_m = rayleigh(1);
beta_k = rayleigh(2);

C_rayleigh = alpha_m * M + beta_k * K;
C = C_rayleigh;
C(1, 1) = C(1, 1) + c_iso;

gamma = 0.5;
beta = 0.25;
a0 = 1 / (beta * dt^2);
a1 = gamma / (beta * dt);
a2 = 1 / (beta * dt);
a3 = 1 / (2 * beta) - 1;
a4 = gamma / beta - 1;
a5 = dt * (gamma / (2 * beta) - 1);

u = zeros(4, n);
v = zeros(4, n);
arel = zeros(4, n);
r = ones(4, 1);
p0 = -M * r * ag(1);
arel(:, 1) = M \ (p0 - C * v(:, 1) - K * u(:, 1));

Keff = K + a0 * M + a1 * C;
for i = 1:n-1
    p_next = -M * r * ag(i+1);
    rhs = p_next ...
        + M * (a0 * u(:, i) + a2 * v(:, i) + a3 * arel(:, i)) ...
        + C * (a1 * u(:, i) + a4 * v(:, i) + a5 * arel(:, i));
    u(:, i+1) = Keff \ rhs;
    arel(:, i+1) = a0 * (u(:, i+1) - u(:, i)) - a2 * v(:, i) - a3 * arel(:, i);
    v(:, i+1) = v(:, i) + dt * ((1 - gamma) * arel(:, i) + gamma * arel(:, i+1));
end

top_displacement = u(4, :)';
isolation_displacement = u(1, :)';
top_relative_acceleration = arel(4, :)';
top_absolute_acceleration = top_relative_acceleration + ag;

writematrix([time, ag], fullfile(script_dir, 'normalized_ground_motion_mps2.txt'), 'Delimiter', 'tab');
writematrix([time, top_displacement], fullfile(script_dir, 'topStoDis.txt'), 'Delimiter', 'tab');
writematrix([time, top_absolute_acceleration], fullfile(script_dir, 'topStoAcc.txt'), 'Delimiter', 'tab');
writematrix([time, isolation_displacement], fullfile(script_dir, 'isolationDisplacement.txt'), 'Delimiter', 'tab');

out = struct();
out.software = 'MATLAB';
out.status = 'success';
out.model = struct();
out.model.units = 'SI';
out.model.masses_kg = masses';
out.model.stiffnesses_N_per_m = stiffnesses';
out.model.isolation_damping_Ns_per_m = c_iso;
out.model.rayleigh_damping_ratio = zeta;
out.model.rayleigh_alpha_mass = alpha_m;
out.model.rayleigh_beta_stiffness = beta_k;
out.model.modal_frequencies_rad_per_s = omega';
out.integration = struct('method', 'Newmark average acceleration', 'gamma', gamma, 'beta', beta);
out.input = struct();
out.input.source = input_path;
out.input.raw_peak = raw_peak;
out.input.target_pga_mps2 = target_pga;
out.input.pga_after_normalization_mps2 = max(abs(ag));
out.input.pga_after_normalization_g = max(abs(ag)) / g;
out.time_step = dt;
out.number_of_samples = n;
out.time = time';
out.ground_acceleration_mps2 = ag';
out.top_story_displacement_m = top_displacement';
out.top_story_acceleration_mps2 = top_absolute_acceleration';
out.top_story_relative_acceleration_mps2 = top_relative_acceleration';
out.isolation_layer_displacement_m = isolation_displacement';

fid = fopen(fullfile(script_dir, 'response.json'), 'w');
fprintf(fid, '%s', jsonencode(out));
fclose(fid);

fprintf('MATLAB isolated-building analysis completed.\n');
fprintf('Samples: %d, dt: %.6f s, normalized PGA: %.12g m/s^2\n', n, dt, max(abs(ag)));
fprintf('Peak top displacement: %.12g m\n', max(abs(top_displacement)));
fprintf('Peak top absolute acceleration: %.12g m/s^2\n', max(abs(top_absolute_acceleration)));
fprintf('Peak isolation displacement: %.12g m\n', max(abs(isolation_displacement)));
