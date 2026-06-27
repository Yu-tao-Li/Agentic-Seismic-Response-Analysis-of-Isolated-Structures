function run_matlab_frame()
%RUN_MATLAB_FRAME Reference 2D isolated frame benchmark for Example 6.

script_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(script_dir));
log_file = fullfile(script_dir, 'execution_log.txt');
if exist(log_file, 'file')
    delete(log_file);
end
diary(log_file);
cleanup_obj = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('MATLAB 2D frame Bouc-Wen benchmark started: %s\n', datestr(now, 31));

p = model_parameters();
input_file = fullfile(root_dir, 'input_data', 'Northridge_01_NO_968.txt');
raw_ag = readmatrix(input_file);
raw_ag = raw_ag(:);
if isempty(raw_ag)
    error('Ground-motion file is empty: %s', input_file);
end
pga_raw = max(abs(raw_ag));
scale = p.target_pga / pga_raw;
ag = raw_ag * scale;
n_steps = numel(ag);
n_intervals = n_steps - 1;
t = (0:(n_steps - 1))' * p.dt;
fprintf('Input samples: %d, raw PGA: %.12g, scale: %.12g, normalized PGA: %.12g m/s^2\n', ...
    n_steps, pga_raw, scale, max(abs(ag)));

[K_frame, K_linear, M, dof] = assemble_reduced_frame(p);
K_initial = K_linear;
for j = 1:3
    K_initial(dof.base_ux(j), dof.base_ux(j)) = K_initial(dof.base_ux(j), dof.base_ux(j)) + p.k0;
end

[periods, omegas] = condensed_modes(K_initial, M, 2);
[rayleigh_alpha_m, rayleigh_beta_k] = rayleigh_coefficients(omegas(1), omegas(2), p.damping_ratio);
C = rayleigh_alpha_m * M + rayleigh_beta_k * K_initial;
fprintf('Elastic modal periods used for damping: %.12g s, %.12g s\n', periods(1), periods(2));
fprintf('Rayleigh coefficients: alphaM = %.12g, betaKinit = %.12g\n', rayleigh_alpha_m, rayleigh_beta_k);

n = size(M, 1);
r = zeros(n, 1);
r(dof.all_ux) = 1.0;

u = zeros(n, 1);
v = zeros(n, 1);
a = zeros(n, 1);
z = zeros(3, 1);

u_hist = zeros(n_steps, n);
roof_hist = zeros(n_steps, 1);
floor_hist = zeros(n_steps, 3);
base_hist = zeros(n_steps, 3);
drift_hist = zeros(n_steps, 3);
iso_disp_hist = zeros(n_steps, 3);
iso_force_hist = zeros(n_steps, 3);
z_hist = zeros(n_steps, 3);
base_shear_hist = zeros(n_steps, 1);
iter_hist = zeros(n_intervals, 1);
residual_hist = zeros(n_intervals, 1);

[iso_force, ~, z_eval] = isolator_forces_and_tangent(u, u, z, dof, p);
iso_force_hist(1, :) = iso_force(:)';
z_hist(1, :) = z_eval(:)';
[roof_hist(1), floor_hist(1, :), base_hist(1, :), drift_hist(1, :), iso_disp_hist(1, :), base_shear_hist(1)] = ...
    response_quantities(u, iso_force, dof, p);

newmark_beta = 0.25;
newmark_gamma = 0.5;
a0 = 1.0 / (newmark_beta * p.dt^2);
a1 = newmark_gamma / (newmark_beta * p.dt);
max_iter = 35;
tol = 1.0e-8;
failed_steps = 0;
max_ratio = 0.0;

for step = 1:n_intervals
    u_old = u;
    v_old = v;
    a_old = a;
    z_old = z;

    u_pred = u_old + p.dt * v_old + p.dt^2 * (0.5 - newmark_beta) * a_old;
    v_const = v_old + p.dt * (1.0 - newmark_gamma) * a_old;
    ext = -M * r * ag(step + 1);

    u_trial = u_old;
    converged = false;
    ratio = Inf;
    for iter = 1:max_iter
        a_trial = a0 * (u_trial - u_pred);
        v_trial = v_const + newmark_gamma * p.dt * a_trial;
        [f_iso, kt_iso, z_trial] = isolator_forces_and_tangent(u_trial, u_old, z_old, dof, p);
        fint = K_linear * u_trial;
        for j = 1:3
            fint(dof.base_ux(j)) = fint(dof.base_ux(j)) + f_iso(j);
        end
        residual = M * a_trial + C * v_trial + fint - ext;
        denom = max([1.0, norm(ext, inf), norm(fint, inf)]);
        ratio = norm(residual, inf) / denom;
        if ratio < tol
            converged = true;
            break;
        end
        K_tangent = K_linear;
        for j = 1:3
            K_tangent(dof.base_ux(j), dof.base_ux(j)) = K_tangent(dof.base_ux(j), dof.base_ux(j)) + kt_iso(j);
        end
        K_eff = K_tangent + a0 * M + a1 * C;
        du = -K_eff \ residual;
        u_trial = u_trial + du;
    end

    if ~converged
        failed_steps = failed_steps + 1;
        warning('Step %d did not converge: residual ratio %.12g after %d iterations', step, ratio, max_iter);
    end

    a = a0 * (u_trial - u_pred);
    v = v_const + newmark_gamma * p.dt * a;
    u = u_trial;
    [iso_force, ~, z] = isolator_forces_and_tangent(u, u_old, z_old, dof, p);

    u_hist(step + 1, :) = u(:)';
    iso_force_hist(step + 1, :) = iso_force(:)';
    z_hist(step + 1, :) = z(:)';
    iter_hist(step) = iter;
    residual_hist(step) = ratio;
    max_ratio = max(max_ratio, ratio);
    [roof_hist(step + 1), floor_hist(step + 1, :), base_hist(step + 1, :), drift_hist(step + 1, :), ...
        iso_disp_hist(step + 1, :), base_shear_hist(step + 1)] = response_quantities(u, iso_force, dof, p);
end

loop_area = zeros(1, 3);
for j = 1:3
    loop_area(j) = sum(0.5 * (iso_force_hist(2:end, j) + iso_force_hist(1:end-1, j)) .* ...
        (iso_disp_hist(2:end, j) - iso_disp_hist(1:end-1, j)));
end

peaks = struct();
peaks.roof_displacement_abs_max = max(abs(roof_hist));
peaks.max_interstory_drift_ratio_abs = max(abs(drift_hist), [], 'all');
peaks.isolation_displacement_abs_max = max(abs(iso_disp_hist), [], 'all');
peaks.mean_isolation_displacement_abs_max = max(abs(mean(iso_disp_hist, 2)));
peaks.total_base_shear_abs_max = max(abs(base_shear_hist));
peaks.isolator_loop_area = loop_area;

result = struct();
result.software = 'MATLAB';
result.status = 'OK';
result.generated_in_this_run = true;
result.time_step = p.dt;
result.num_input_samples = n_steps;
result.num_analysis_intervals = n_intervals;
result.input_file = relative_path(input_file, root_dir);
result.raw_pga = pga_raw;
result.normalization_scale = scale;
result.normalized_pga = max(abs(ag));
result.fundamental_period = periods(1);
result.elastic_periods = periods(:)';
result.damping = struct('ratio_modes_1_2', p.damping_ratio, ...
    'rayleigh_alpha_m', rayleigh_alpha_m, 'rayleigh_beta_k_initial', rayleigh_beta_k, ...
    'stiffness_matrix', 'initial elastic tangent including k0 horizontal isolator tangent');
result.bouc_wen_update = ['Backward-Euler increment equation: z_{n+1}-z_n-Ao*du+', ...
    'beta*abs(du)*abs(z_{n+1})^(n-1)*z_{n+1}+gamma*du*abs(z_{n+1})^n = 0; ', ...
    'F = alpha*k0*u + (1-alpha)*k0*z_{n+1}.'];
result.modeling_assumptions = assumptions_struct();
result.time = t;
result.ground_acceleration = ag;
result.roof_displacement = roof_hist;
result.floor_displacements = floor_hist;
result.base_node_displacements = base_hist;
result.interstory_drift_ratios = drift_hist;
result.isolation_displacement = mean(iso_disp_hist, 2);
result.isolator_displacements = iso_disp_hist;
result.total_base_shear = base_shear_hist;
result.isolator_forces = iso_force_hist;
result.bouc_wen_z = z_hist;
result.hysteresis_loop_area = loop_area;
result.peak_responses = peaks;
result.convergence = struct('relative_residual_tolerance', tol, 'max_iterations_per_step', max_iter, ...
    'failed_step_count', failed_steps, 'max_relative_residual', max_ratio, ...
    'max_iterations_observed', max(iter_hist), 'iteration_history', iter_hist, ...
    'relative_residual_history', residual_hist);

json_file = fullfile(script_dir, 'matlab_results.json');
fid = fopen(json_file, 'w');
fprintf(fid, '%s', jsonencode(result));
fclose(fid);

time_table = table(t, ag, roof_hist, drift_hist(:, 1), drift_hist(:, 2), drift_hist(:, 3), ...
    mean(iso_disp_hist, 2), base_shear_hist, ...
    'VariableNames', {'time_s', 'ground_accel_mps2', 'roof_disp_m', 'drift_story1', 'drift_story2', ...
    'drift_story3', 'mean_isolation_disp_m', 'total_base_shear_N'});
writetable(time_table, fullfile(script_dir, 'time_histories.csv'));

hyst_table = table(t, iso_disp_hist(:, 1), iso_force_hist(:, 1), z_hist(:, 1), ...
    iso_disp_hist(:, 2), iso_force_hist(:, 2), z_hist(:, 2), ...
    iso_disp_hist(:, 3), iso_force_hist(:, 3), z_hist(:, 3), ...
    'VariableNames', {'time_s', 'u_iso_1_m', 'f_iso_1_N', 'z_iso_1_m', ...
    'u_iso_2_m', 'f_iso_2_N', 'z_iso_2_m', 'u_iso_3_m', 'f_iso_3_N', 'z_iso_3_m'});
writetable(hyst_table, fullfile(script_dir, 'hysteresis.csv'));

fprintf('MATLAB benchmark completed: failed_steps=%d, max residual ratio=%.12g\n', failed_steps, max_ratio);
fprintf('Outputs written to %s\n', script_dir);
end

function p = model_parameters()
p.E = 2.06e11;
p.nu = 0.30;
p.rho = 7850.0;
p.Ac = 0.020;
p.Ic = 8.0e-4;
p.Ab = 0.015;
p.Ib = 4.5e-4;
p.bay_width = 6.0;
p.story_height = 3.6;
p.floor_mass = 2.5e5;
p.k0 = 2.0e6;
p.Fy = 1.0e4;
p.uy = p.Fy / p.k0;
p.alpha = 0.05;
p.bouc_n = 2.0;
p.beta_bw = 2.0e4;
p.gamma_bw = 2.0e4;
p.Ao = 1.0;
p.deltaA = 0.0;
p.deltaNu = 0.0;
p.deltaEta = 0.0;
p.kv = 1.0e10;
p.g = 9.81;
p.target_pga = 0.40 * p.g;
p.dt = 0.01;
p.damping_ratio = 0.05;
end

function [K_red_frame, K_red_linear, M_red, dof] = assemble_reduced_frame(p)
n_cols = 3;
n_levels = 4;
n_nodes = n_cols * n_levels;
n_dof_full = n_nodes * 3;
coords = zeros(n_nodes, 2);
for lev = 0:(n_levels - 1)
    for col = 0:(n_cols - 1)
        node = lev * n_cols + col + 1;
        coords(node, :) = [col * p.bay_width, lev * p.story_height];
    end
end

K_full = zeros(n_dof_full);
for col = 0:(n_cols - 1)
    for lev = 0:2
        ni = lev * n_cols + col + 1;
        nj = (lev + 1) * n_cols + col + 1;
        K_full = add_frame_element(K_full, coords, ni, nj, p.E, p.Ac, p.Ic);
    end
end
for lev = 1:3
    for col = 0:1
        ni = lev * n_cols + col + 1;
        nj = lev * n_cols + col + 2;
        K_full = add_frame_element(K_full, coords, ni, nj, p.E, p.Ab, p.Ib);
    end
end

M_full = zeros(n_dof_full);
for lev = 1:3
    for col = 0:2
        node = lev * n_cols + col + 1;
        M_full(node_dofs(node, 1), node_dofs(node, 1)) = p.floor_mass / 3.0;
    end
end

map = zeros(n_dof_full, 1);
count = 0;
floor_master = zeros(3, 1);
for lev = 0:3
    if lev >= 1
        count = count + 1;
        floor_master(lev) = count;
        for col = 0:2
            node = lev * n_cols + col + 1;
            map(node_dofs(node, 1)) = floor_master(lev);
        end
    else
        for col = 0:2
            node = col + 1;
            count = count + 1;
            map(node_dofs(node, 1)) = count;
        end
    end
    for col = 0:2
        node = lev * n_cols + col + 1;
        count = count + 1;
        map(node_dofs(node, 2)) = count;
        count = count + 1;
        map(node_dofs(node, 3)) = count;
    end
end

T = zeros(n_dof_full, count);
for i = 1:n_dof_full
    T(i, map(i)) = 1.0;
end

K_red_frame = T' * K_full * T;
M_red = T' * M_full * T;
K_red_linear = K_red_frame;

base_ux = zeros(3, 1);
base_uy = zeros(3, 1);
for col = 0:2
    node = col + 1;
    base_ux(col + 1) = map(node_dofs(node, 1));
    base_uy(col + 1) = map(node_dofs(node, 2));
    K_red_linear(base_uy(col + 1), base_uy(col + 1)) = K_red_linear(base_uy(col + 1), base_uy(col + 1)) + p.kv;
end

floor_ux = floor_master(:);
all_ux = unique(map(1:3:end), 'stable');

dof = struct();
dof.map = map;
dof.floor_ux = floor_ux;
dof.base_ux = base_ux;
dof.base_uy = base_uy;
dof.roof_ux = floor_ux(3);
dof.all_ux = all_ux;
end

function K = add_frame_element(K, coords, ni, nj, E, A, I)
xi = coords(ni, 1);
yi = coords(ni, 2);
xj = coords(nj, 1);
yj = coords(nj, 2);
dx = xj - xi;
dy = yj - yi;
L = sqrt(dx^2 + dy^2);
c = dx / L;
s = dy / L;
EA_L = E * A / L;
EI = E * I;
k = [ EA_L,       0,            0, -EA_L,       0,            0;
          0, 12*EI/L^3,  6*EI/L^2,      0, -12*EI/L^3,  6*EI/L^2;
          0,  6*EI/L^2,     4*EI/L,      0,  -6*EI/L^2,     2*EI/L;
      -EA_L,       0,            0,  EA_L,       0,            0;
          0,-12*EI/L^3, -6*EI/L^2,      0,  12*EI/L^3, -6*EI/L^2;
          0,  6*EI/L^2,     2*EI/L,      0,  -6*EI/L^2,     4*EI/L];
R = [ c, s, 0, 0, 0, 0;
     -s, c, 0, 0, 0, 0;
      0, 0, 1, 0, 0, 0;
      0, 0, 0, c, s, 0;
      0, 0, 0,-s, c, 0;
      0, 0, 0, 0, 0, 1];
kg = R' * k * R;
dofs = [node_dofs(ni, 1), node_dofs(ni, 2), node_dofs(ni, 3), ...
        node_dofs(nj, 1), node_dofs(nj, 2), node_dofs(nj, 3)];
K(dofs, dofs) = K(dofs, dofs) + kg;
end

function d = node_dofs(node, comp)
d = (node - 1) * 3 + comp;
end

function [periods, omegas] = condensed_modes(K, M, n_modes)
mass_diag = diag(M);
m_dofs = find(mass_diag > 0);
z_dofs = find(mass_diag == 0);
Kmm = K(m_dofs, m_dofs);
Kmz = K(m_dofs, z_dofs);
Kzz = K(z_dofs, z_dofs);
Kcond = Kmm - Kmz * (Kzz \ Kmz');
Mcond = M(m_dofs, m_dofs);
[V, D] = eig((Kcond + Kcond') / 2, (Mcond + Mcond') / 2); %#ok<ASGLU>
lambda = real(diag(D));
lambda = lambda(lambda > 0);
omega_all = sort(sqrt(lambda));
omegas = omega_all(1:n_modes);
periods = 2 * pi ./ omegas;
end

function [alpha_m, beta_k] = rayleigh_coefficients(w1, w2, xi)
A = [1.0 / (2.0 * w1), w1 / 2.0; 1.0 / (2.0 * w2), w2 / 2.0];
x = A \ [xi; xi];
alpha_m = x(1);
beta_k = x(2);
end

function [force, tangent, z_new] = isolator_forces_and_tangent(u_trial, u_old, z_old, dof, p)
force = zeros(3, 1);
tangent = zeros(3, 1);
z_new = zeros(3, 1);
for j = 1:3
    idx = dof.base_ux(j);
    du = u_trial(idx) - u_old(idx);
    [zj, dzdu] = boucwen_backward_euler(z_old(j), du, p);
    uj = u_trial(idx);
    force(j) = p.alpha * p.k0 * uj + (1.0 - p.alpha) * p.k0 * zj;
    tangent(j) = p.alpha * p.k0 + (1.0 - p.alpha) * p.k0 * dzdu;
    z_new(j) = zj;
end
end

function [z, dzdu] = boucwen_backward_euler(z_old, du, p)
z = z_old + p.Ao * du;
for it = 1:30
    absz = abs(z);
    if absz < 1.0e-14
        absz = 0.0;
    end
    g = z - z_old - p.Ao * du + p.beta_bw * abs(du) * abs(z)^(p.bouc_n - 1.0) * z + ...
        p.gamma_bw * du * abs(z)^p.bouc_n;
    if abs(g) < 1.0e-13 * max(1.0, abs(z))
        break;
    end
    if absz == 0.0
        dgdz = 1.0;
    else
        dgdz = 1.0 + p.beta_bw * abs(du) * p.bouc_n * absz^(p.bouc_n - 1.0) + ...
            p.gamma_bw * du * p.bouc_n * absz^(p.bouc_n - 1.0) * sign(z);
    end
    z = z - g / dgdz;
end
absz = abs(z);
if absz == 0.0
    dgdz = 1.0;
else
    dgdz = 1.0 + p.beta_bw * abs(du) * p.bouc_n * absz^(p.bouc_n - 1.0) + ...
        p.gamma_bw * du * p.bouc_n * absz^(p.bouc_n - 1.0) * sign(z);
end
dgdu = -p.Ao + p.beta_bw * sign_nonzero(du) * abs(z)^(p.bouc_n - 1.0) * z + ...
    p.gamma_bw * abs(z)^p.bouc_n;
dzdu = -dgdu / dgdz;
end

function s = sign_nonzero(x)
if x > 0
    s = 1.0;
elseif x < 0
    s = -1.0;
else
    s = 0.0;
end
end

function [roof, floors, bases, drifts, iso_disp, base_shear] = response_quantities(u, iso_force, dof, p)
floors = [u(dof.floor_ux(1)), u(dof.floor_ux(2)), u(dof.floor_ux(3))];
bases = [u(dof.base_ux(1)), u(dof.base_ux(2)), u(dof.base_ux(3))];
base_avg = mean(bases);
roof = floors(3);
drifts = [(floors(1) - base_avg) / p.story_height, ...
          (floors(2) - floors(1)) / p.story_height, ...
          (floors(3) - floors(2)) / p.story_height];
iso_disp = bases;
base_shear = sum(iso_force);
end

function s = assumptions_struct()
s = struct();
s.frame = '2D Euler-Bernoulli elastic beam-column stiffness, 3 DOF per node before diaphragm reduction.';
s.mass = 'Only specified lumped floor mass is assigned to horizontal floor DOFs; steel density is recorded but not added as self-mass.';
s.constraints = 'At each floor, the three horizontal DOFs are transformed to one rigid-diaphragm coordinate; vertical translations and rotations remain independent.';
s.isolation = 'Three horizontal Bouc-Wen isolators and three vertical linear springs connect top isolator nodes to fixed ground.';
s.story_drift = 'Story 1 drift uses floor-1 diaphragm displacement minus average top-of-isolator displacement divided by story height.';
s.base_shear = 'Total base shear is the sum of the three positive-resisting horizontal isolator forces.';
s.coordinates = 'All output displacements are relative to ground in the frame x direction.';
end

function rp = relative_path(pathname, root_dir)
rp = char(java.io.File(pathname).toURI().relativize(java.io.File(pathname).toURI()));
if isempty(rp)
    rp = strrep(pathname, [root_dir filesep], '');
end
end
