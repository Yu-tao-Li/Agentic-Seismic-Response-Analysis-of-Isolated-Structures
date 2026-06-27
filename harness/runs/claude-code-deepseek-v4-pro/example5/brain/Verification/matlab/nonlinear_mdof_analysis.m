% Nonlinear MDOF Seismic Analysis - Example 5
% 4-story seismically isolated structure
% Isolation layer (story 1): bilinear kinematic hardening (b=0.05)
% Upper stories (2-4): ideal elastic-plastic (b=0)
% Newmark-beta integration with Newton-Raphson iteration

clear; clc;

% ===== Load seismic data =====
data = load('../../input_data/Northridge_01_NO_968.txt');
dt = 0.01;
g = 9.81;
PGA = 0.40 * g;

acceleration_data = data(:,1);
time = (0:length(acceleration_data)-1)' * dt;
acceleration_data = acceleration_data / max(abs(acceleration_data)) * PGA;

% ===== Structure parameters =====
m = [250 270 270 180] * 1e3;      % kg
k = [50 245 195 98] * 1e6;         % N/m
fy = [500 1225 975 490] * 1e3;     % N
b = [0.05 0 0 0];                  % post-yield hardening ratio
n_stories = 4;

% ===== Rayleigh damping =====
damping_ratio = 0.05;
M = diag(m);
K0 = diag([k(1)+k(2), k(2)+k(3), k(3)+k(4), k(4)]) ...
     - diag(k(2:end), 1) - diag(k(2:end), -1);

[V, D] = eig(K0, M);
omega = sqrt(diag(D));
alpha = 2 * damping_ratio * omega(1) * omega(2) / (omega(1) + omega(2));
beta = 2 * damping_ratio / (omega(1) + omega(2));
C = alpha * M + beta * K0;

% ===== Newmark-beta parameters =====
beta_n = 0.25;
gamma = 0.5;

% ===== Initial conditions =====
u = zeros(n_stories, 1);
u_dot = zeros(n_stories, 1);
Rs = zeros(n_stories, 1);
I_vec = ones(n_stories, 1);
u_dot_dot = M \ (-acceleration_data(1) * M * I_vec - C * u_dot - Rs);

% ===== Storage arrays =====
n_steps = length(time);
top_displacement = zeros(n_steps, 1);
top_acceleration = zeros(n_steps, 1);
iso_displacement = zeros(n_steps, 1);
base_shear = zeros(n_steps, 1);
hyst_force = zeros(n_steps, 1);
hyst_disp = zeros(n_steps, 1);
yield_flags = zeros(n_steps, n_stories);
conv_fails = zeros(n_steps, 1);
iter_count = zeros(n_steps, 1);

% Initial step storage
top_displacement(1) = u(end);
top_acceleration(1) = u_dot_dot(end);
iso_displacement(1) = u(1);
base_shear(1) = sum(Rs);

% Initial effective stiffness (using elastic K)
K_eff = K0 + (1/(beta_n*dt^2))*M + (gamma/(beta_n*dt))*C;
k_story = k;  % tangent story stiffnesses
k_trial = k;

max_iter = 100;
tol = 1e-5;

fprintf('Starting nonlinear MDOF analysis: %d steps\n', n_steps);

for j = 1:n_steps-1
    % External load increment
    delta_P = -(acceleration_data(j+1) - acceleration_data(j)) * M * I_vec;

    % Effective load increment
    delta_P_eff = delta_P ...
        + ((1/(beta_n*dt))*M + (gamma/beta_n)*C) * u_dot ...
        + ((1/(2*beta_n))*M + dt*(gamma/(2*beta_n) - 1)*C) * u_dot_dot;

    % Newton-Raphson iteration
    u_trial = u;
    Rs_trial = Rs;
    delta_R_trial = delta_P_eff;

    converged = false;
    for iter = 1:max_iter
        if sum(abs(delta_R_trial)) < tol
            converged = true;
            break;
        end

        delta_u = K_eff \ delta_R_trial;
        u_trial_n = u_trial + delta_u;

        % Compute restoring force and tangent stiffness
        [Rs_trial_n, k_trial] = computeRsStiffness(Rs_trial, u_trial, u_trial_n, fy, k, b);
        K_trial = computeStiffnessMatrix(k_trial);
        K_eff_trial = K_trial + (1/(beta_n*dt^2))*M + (gamma/(beta_n*dt))*C;

        delta_F = Rs_trial_n - Rs_trial ...
            + (gamma/(beta_n*dt))*C*delta_u ...
            + (1/(beta_n*dt^2))*M*delta_u;
        delta_R_trial = delta_R_trial - delta_F;

        u_trial = u_trial_n;
        Rs_trial = Rs_trial_n;
        K_eff = K_eff_trial;
    end

    iter_count(j+1) = iter;
    conv_fails(j+1) = ~converged;

    % Update state
    delta_u = u_trial - u;
    delta_u_dot = (gamma/(beta_n*dt))*delta_u - (gamma/beta_n)*u_dot ...
        + (1 - gamma/(2*beta_n))*dt*u_dot_dot;
    delta_u_dot_dot = (1/(beta_n*dt^2))*delta_u - (1/(beta_n*dt))*u_dot ...
        - (1/(2*beta_n))*u_dot_dot;

    u = u_trial;
    u_dot = u_dot + delta_u_dot;
    u_dot_dot = u_dot_dot + delta_u_dot_dot;
    Rs = Rs_trial;
    k_story = k_trial;

    % Store results
    top_displacement(j+1) = u(end);
    top_acceleration(j+1) = u_dot_dot(end);
    iso_displacement(j+1) = u(1);

    % Base shear = sum of all restoring forces = shear at story 1
    base_shear(j+1) = sum(Rs);

    % Hysteresis: isolation layer force vs displacement
    hyst_force(j+1) = sum(Rs);       % total base shear = isolation shear
    hyst_disp(j+1) = u(1);           % isolation displacement

    % Yielding indicators: k_story(i) < k(i) means yielded
    for iy = 1:n_stories
        if k_story(iy) < k(iy) * 0.99
            yield_flags(j+1, iy) = 1;
        end
    end
end

fprintf('Analysis complete.\n');
fprintf('Peak top displacement: %.6f m\n', max(abs(top_displacement)));
fprintf('Peak top acceleration: %.4f m/s^2\n', max(abs(top_acceleration)));
fprintf('Peak isolation displacement: %.6f m\n', max(abs(iso_displacement)));
fprintf('Peak base shear: %.2f N\n', max(abs(base_shear)));
fprintf('Convergence failures: %d/%d\n', sum(conv_fails), n_steps);

% ===== Save output files =====
% topStoDisIso2.txt
topStoDisIso2 = [time, top_displacement];
save('topStoDisIso2.txt', 'topStoDisIso2', '-ascii', '-double');

% topStoAccIso2.txt
topStoAccIso2 = [time, top_acceleration];
save('topStoAccIso2.txt', 'topStoAccIso2', '-ascii', '-double');

% Isolation displacement
iso_disp_out = [time, iso_displacement];
save('iso_displacement.txt', 'iso_disp_out', '-ascii', '-double');

% Base shear
base_shear_out = [time, base_shear];
save('base_shear.txt', 'base_shear_out', '-ascii', '-double');

% Hysteresis data
hysteresis_out = [time, hyst_disp, hyst_force];
save('hysteresis_layer1.txt', 'hysteresis_out', '-ascii', '-double');

% Yielding indicators
yield_out = [time, yield_flags];
save('yield_flags.txt', 'yield_out', '-ascii', '-double');

% Convergence log
conv_out = [time, conv_fails, iter_count];
save('convergence_log.txt', 'conv_out', '-ascii', '-double');

% MATLAB response metadata
response = struct();
response.software = 'MATLAB';
response.version = version();
response.model_type = 'nonlinear_MDOF_Newmark_NR';
response.n_stories = n_stories;
response.parameters = struct('mass', m, 'stiffness', k, 'yield_force', fy, 'hardening_ratio', b);
response.damping = struct('ratio', damping_ratio, 'alpha', alpha, 'beta', beta);
response.integration = struct('method', 'Newmark-beta', 'beta', beta_n, 'gamma', gamma);
response.input = struct('file', 'Northridge_01_NO_968.txt', 'n_rows', length(acceleration_data), 'dt', dt, 'PGA', PGA);
response.results = struct(...
    'peak_top_displacement_m', max(abs(top_displacement)), ...
    'peak_top_acceleration_ms2', max(abs(top_acceleration)), ...
    'peak_isolation_displacement_m', max(abs(iso_displacement)), ...
    'peak_base_shear_N', max(abs(base_shear)), ...
    'convergence_failures', sum(conv_fails), ...
    'convergence_rate_pct', (1-sum(conv_fails)/n_steps)*100);
fid = fopen('matlab_response.json', 'w');
fprintf(fid, '%s', jsonencode(response, 'PrettyPrint', true));
fclose(fid);

fprintf('Output files saved.\n');

% ===== Helper functions =====
function [Rs, stiffness] = computeRsStiffness(lastRs, lastDisp, disp, fy, k0, b)
    n = length(lastRs);
    lastShear = zeros(n, 1);
    lastRelativeDisp = zeros(n, 1);
    relativeDisp = zeros(n, 1);
    shear = zeros(n, 1);
    stiffness = zeros(n, 1);

    for i = 1:n
        lastShear(i) = sum(lastRs(i:n));
    end

    lastRelativeDisp(1) = lastDisp(1);
    relativeDisp(1) = disp(1);
    for i = 2:n
        lastRelativeDisp(i) = lastDisp(i) - lastDisp(i-1);
        relativeDisp(i) = disp(i) - disp(i-1);
    end

    for i = 1:n
        fy_minus_b = fy(i) * (1 - b(i));
        k_sh = b(i) * k0(i);
        deltaDisp = relativeDisp(i) - lastRelativeDisp(i);
        c = lastShear(i) + k0(i) * deltaDisp;
        shear(i) = max(k_sh * relativeDisp(i) - fy_minus_b, ...
                       min(k_sh * relativeDisp(i) + fy_minus_b, c));
        if abs(shear(i) - c) < 1e-3
            stiffness(i) = k0(i);
        else
            stiffness(i) = k_sh;
        end
    end

    Rs = zeros(n, 1);
    Rs(n) = shear(n);
    for i = 1:n-1
        Rs(i) = shear(i) - shear(i+1);
    end
end

function K = computeStiffnessMatrix(storyStiffness)
    n = length(storyStiffness);
    K = zeros(n, n);
    for i = 1:n
        if i < n
            K(i, i) = storyStiffness(i) + storyStiffness(i + 1);
            K(i, i + 1) = -storyStiffness(i + 1);
            K(i + 1, i) = -storyStiffness(i + 1);
        else
            K(i, i) = storyStiffness(i);
        end
    end
end
