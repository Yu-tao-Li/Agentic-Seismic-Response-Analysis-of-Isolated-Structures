%% 2D Three-Story Two-Bay Base-Isolated Steel Frame
%% Bouc-Wen Isolators | Newmark-beta Integration | Newton-Raphson
%  Harness: example6 | Software: MATLAB

clear; close all;

%% ========================================================================
%  1. MATERIAL AND SECTION PROPERTIES
%  ========================================================================
E   = 2.06e11;      % Elastic modulus [Pa]
nu  = 0.30;         % Poisson ratio
rho = 7850;         % Density [kg/m^3]
A_c = 0.020;       % Column area [m^2]
I_c = 8.0e-4;      % Column 2nd moment of area [m^4]
A_b = 0.015;       % Beam area [m^2]
I_b = 4.5e-4;      % Beam 2nd moment of area [m^4]

%% ========================================================================
%  2. GEOMETRY AND NODE DEFINITION
%  ========================================================================
n_stories = 3;
n_bays    = 2;
H_story   = 3.6;    % Story height [m]
L_bay     = 6.0;    % Bay width [m]

% Node numbering (12 superstructure nodes):
%   Level 0 (isolation top): nodes 1,2,3
%   Level 1 (floor 1):       nodes 4,5,6
%   Level 2 (floor 2):       nodes 7,8,9
%   Level 3 (roof):          nodes 10,11,12

n_nodes    = 12;
ndf        = 3;     % DOFs per node: ux, uy, rz
n_dofs    = n_nodes * ndf;  % 36

% Node coordinates [x, y]
coords = zeros(n_nodes, 2);
node_id = 1;
for level = 0:n_stories
    y = level * H_story;
    for col = 1:(n_bays+1)
        x = (col-1) * L_bay;
        coords(node_id, :) = [x, y];
        node_id = node_id + 1;
    end
end

%% ========================================================================
%  3. ELEMENT DEFINITION
%  ========================================================================
% Build elements
elements = struct('type',{},'nodes',{},'L',{},'A',{},'I',{},'Ke_global',{});

% Columns: 3 columns per story
elem_count = 0;
for story = 0:(n_stories-1)
    for col = 1:(n_bays+1)
        n1 = story*(n_bays+1) + col;
        n2 = (story+1)*(n_bays+1) + col;
        elem_count = elem_count + 1;
        elements(elem_count).type = 'column';
        elements(elem_count).nodes = [n1, n2];
        elements(elem_count).A = A_c;
        elements(elem_count).I = I_c;
    end
end

% Beams: 2 beams per floor level
for level = 1:n_stories
    for bay = 1:n_bays
        n1 = level*(n_bays+1) + bay;
        n2 = level*(n_bays+1) + bay + 1;
        elem_count = elem_count + 1;
        elements(elem_count).type = 'beam';
        elements(elem_count).nodes = [n1, n2];
        elements(elem_count).A = A_b;
        elements(elem_count).I = I_b;
    end
end

n_elements = elem_count;  % 9 columns + 6 beams = 15

%% ========================================================================
%  4. ELEMENT STIFFNESS MATRICES
%  ========================================================================
for e = 1:n_elements
    n1 = elements(e).nodes(1);
    n2 = elements(e).nodes(2);
    dx = coords(n2,1) - coords(n1,1);
    dy = coords(n2,2) - coords(n1,2);
    L = sqrt(dx^2 + dy^2);
    cos_a = dx/L;
    sin_a = dy/L;

    % Local stiffness: DOF order ux1, uy1, rz1, ux2, uy2, rz2
    Ke_local = zeros(6,6);
    EA_L = E * elements(e).A / L;
    EI_L3 = E * elements(e).I / L^3;
    EI_L2 = E * elements(e).I / L^2;
    EI_L  = E * elements(e).I / L;

    Ke_local(1,1) = EA_L;  Ke_local(1,4) = -EA_L;
    Ke_local(4,1) = -EA_L; Ke_local(4,4) = EA_L;
    Ke_local(2,2) = 12*EI_L3;   Ke_local(2,3) = 6*EI_L2;
    Ke_local(2,5) = -12*EI_L3;  Ke_local(2,6) = 6*EI_L2;
    Ke_local(3,2) = 6*EI_L2;    Ke_local(3,3) = 4*EI_L;
    Ke_local(3,5) = -6*EI_L2;   Ke_local(3,6) = 2*EI_L;
    Ke_local(5,2) = -12*EI_L3;  Ke_local(5,3) = -6*EI_L2;
    Ke_local(5,5) = 12*EI_L3;   Ke_local(5,6) = -6*EI_L2;
    Ke_local(6,2) = 6*EI_L2;    Ke_local(6,3) = 2*EI_L;
    Ke_local(6,5) = -6*EI_L2;   Ke_local(6,6) = 4*EI_L;

    % Rotation matrix
    lam = [cos_a, sin_a, 0; -sin_a, cos_a, 0; 0, 0, 1];
    T = [lam, zeros(3); zeros(3), lam];

    elements(e).Ke_global = T' * Ke_local * T;
    elements(e).L = L;
    elements(e).cos_a = cos_a;
    elements(e).sin_a = sin_a;
end

%% ========================================================================
%  5. GLOBAL STIFFNESS MATRIX ASSEMBLY
%  ========================================================================
K_global = zeros(n_dofs, n_dofs);
for e = 1:n_elements
    n1 = elements(e).nodes(1);
    n2 = elements(e).nodes(2);
    dofs = [3*(n1-1)+1:3*n1, 3*(n2-1)+1:3*n2];
    K_global(dofs, dofs) = K_global(dofs, dofs) + elements(e).Ke_global;
end

% Save physical frame stiffness BEFORE penalty constraints
K_frame_only = K_global;

%% ========================================================================
%  6. MASS MATRIX (lumped, horizontal DOF only)
%  ========================================================================
M_floor = 2.5e5;  % Total floor mass [kg]
m_node = M_floor / (n_bays+1);  % Mass per floor node

M_global = zeros(n_dofs, n_dofs);
for level = 1:n_stories
    for col = 1:(n_bays+1)
        node = level*(n_bays+1) + col;
        dof_ux = 3*(node-1) + 1; % horizontal DOF only
        M_global(dof_ux, dof_ux) = M_global(dof_ux, dof_ux) + m_node;
    end
end

%% ========================================================================
%  7. RIGID DIAPHRAGM CONSTRAINTS (Penalty Method)
%  ========================================================================
% Add large penalty stiffness to tie horizontal DOFs at each floor level
K_penalty_diaphragm = 1.0e12;  % N/m, >> beam axial stiffness ~5e8

for level = 1:n_stories
    nodes_at_level = level*(n_bays+1) + (1:(n_bays+1));
    % Tie each node to the first node at this level
    master_node = nodes_at_level(1);
    dof_master = 3*(master_node-1) + 1;
    for j = 2:length(nodes_at_level)
        slave_node = nodes_at_level(j);
        dof_slave = 3*(slave_node-1) + 1;
        % Add penalty between master and slave
        K_global(dof_master, dof_master) = K_global(dof_master, dof_master) + K_penalty_diaphragm;
        K_global(dof_master, dof_slave) = K_global(dof_master, dof_slave) - K_penalty_diaphragm;
        K_global(dof_slave, dof_master) = K_global(dof_slave, dof_master) - K_penalty_diaphragm;
        K_global(dof_slave, dof_slave) = K_global(dof_slave, dof_slave) + K_penalty_diaphragm;
    end
end

%% ========================================================================
%  8. ISOLATOR PARAMETERS (Bouc-Wen)
%  ========================================================================
k0    = 2.0e6;     % Initial horizontal stiffness [N/m]
Fy    = 1.0e4;     % Yield force [N]
uy_val = Fy / k0;  % Yield deformation [m]
alpha = 0.05;      % Post-yield stiffness ratio
n_bw  = 2.0;       % Smoothness exponent
beta  = 2.0e4;     % Bouc-Wen shape parameter
gamma_bw = 2.0e4;  % Bouc-Wen shape parameter (rename to avoid conflict)
Ao    = 1.0;       % Reference amplitude
kv    = 1.0e10;    % Vertical stiffness per isolator [N/m]

n_isolators = n_bays + 1;  % 3 isolators

%% ========================================================================
%  9. BUILD CLEAN TANGENT STIFFNESS FOR EIGENVALUE AND DAMPING
%  ========================================================================
% Use physical stiffness only (frame + isolators, NO penalty constraints)
K_tangent = K_frame_only;

% Add isolator elastic stiffness at level 0 nodes
for i_iso = 1:n_isolators
    node = i_iso;  % Level 0 node
    dof_ux = 3*(node-1) + 1;
    K_tangent(dof_ux, dof_ux) = K_tangent(dof_ux, dof_ux) + k0;
    dof_uy = 3*(node-1) + 2;
    K_tangent(dof_uy, dof_uy) = K_tangent(dof_uy, dof_uy) + kv;
end

%% ========================================================================
% 10. MODAL ANALYSIS AND RAYLEIGH DAMPING
%  ========================================================================
[V_k, D_k] = eig(K_tangent, M_global);
omega_sq = real(diag(D_k));
% Filter valid modes
valid = omega_sq > 1e-3;
omega_sq = omega_sq(valid);
omega = sqrt(omega_sq);
omega = sort(omega);

T_periods = 2*pi ./ omega;
freqs = omega / (2*pi);

fprintf('=== Modal Properties (Elastic Initial State) ===\n');
for i = 1:min(6, length(T_periods))
    fprintf('Mode %d: T = %.4f s, f = %.4f Hz\n', i, T_periods(i), freqs(i));
end

% Rayleigh damping: 5% in first two modes
zeta = 0.05;
omega1 = omega(1);
omega2 = omega(2);

A_ray = [1/omega1, omega1; 1/omega2, omega2];
coeffs = A_ray \ [2*zeta; 2*zeta];
alpha_m = coeffs(1);
beta_k  = coeffs(2);

fprintf('\nRayleigh damping: alpha_m=%.6e, beta_k=%.6e\n', alpha_m, beta_k);
fprintf('Target: zeta=%.2f%% at T1=%.4fs, T2=%.4fs\n', zeta*100, T_periods(1), T_periods(2));

% Verify
zet1 = 0.5 * (alpha_m/omega1 + beta_k*omega1);
zet2 = 0.5 * (alpha_m/omega2 + beta_k*omega2);
fprintf('Check: zeta1=%.4f%%, zeta2=%.4f%%\n', zet1*100, zet2*100);

% Damping matrix: use physical stiffness only (no penalty)
C_damp = alpha_m * M_global + beta_k * K_tangent;

%% ========================================================================
% 11. GROUND MOTION INPUT
%  ========================================================================
gm_data_raw = load('harness\runs\claude-code-deepseek-v4-pro\example6\brain\input_data\Northridge_01_NO_968.txt');

g_val = 9.81;
target_pga = 0.40 * g_val;
raw_pga = max(abs(gm_data_raw));
scale_factor = target_pga / raw_pga;
ag_history = gm_data_raw * scale_factor;

fprintf('\nGround motion: raw PGA=%.4f m/s^2, scaled PGA=%.4f m/s^2, scale=%.6f\n', ...
    raw_pga, target_pga, scale_factor);

dt = 0.01;
n_steps = length(ag_history);
t_vec = (0:n_steps-1)' * dt;

%% ========================================================================
% 12. NEWMARK INTEGRATION CONSTANTS
%  ========================================================================
beta_nm = 0.25;
gamma_nm = 0.5;

a0 = 1 / (beta_nm * dt^2);
a1 = gamma_nm / (beta_nm * dt);
a2 = 1 / (beta_nm * dt);
a3 = 0.5 / beta_nm - 1;
a4 = gamma_nm / beta_nm - 1;
a5 = 0.5 * dt * (gamma_nm / beta_nm - 2);

%% ========================================================================
% 13. INITIALIZATION
%  ========================================================================
U = zeros(n_dofs, 1);
V = zeros(n_dofs, 1);
A_acc = zeros(n_dofs, 1);
z_bw = zeros(n_isolators, 1);

% Effective stiffness matrix (constant linear part)
K_eff_linear = K_global + a0 * M_global + a1 * C_damp;

%% ========================================================================
% 14. STORAGE
%  ========================================================================
store_t     = zeros(n_steps, 1);
store_roof_ux    = zeros(n_steps, 1);
store_iso_ux     = zeros(n_steps, 3);
store_iso_F      = zeros(n_steps, 3);
store_iso_z      = zeros(n_steps, 3);
store_base_shear = zeros(n_steps, 1);
store_drift1     = zeros(n_steps, 1);
store_drift2     = zeros(n_steps, 1);
store_drift3     = zeros(n_steps, 1);
store_n_iter     = zeros(n_steps, 1);
store_converged  = zeros(n_steps, 1);

peak_roof_ux = 0; peak_iso_ux = 0; peak_base_shear = 0;
peak_drift1 = 0; peak_drift2 = 0; peak_drift3 = 0;
failed_steps = 0;
tol_nr     = 1e-8;
max_iter   = 50;

%% ========================================================================
% 15. TIME STEPPING
%  ========================================================================
fprintf('\n=== Time Integration (%d steps) ===\n', n_steps);

for step = 1:n_steps
    ag = ag_history(step);

    % Store previous converged displacement
    U_prev = U;

    % Newmark predictor
    U_pred = U + dt * V + 0.5*dt^2*(1-2*beta_nm)*A_acc;
    V_pred = V + dt * (1-gamma_nm)*A_acc;

    % Initial guess for incremental displacement
    delta_U = zeros(n_dofs, 1);
    A_new = zeros(n_dofs, 1);

    converged = false;
    n_iter = 0;

    while ~converged && n_iter < max_iter
        n_iter = n_iter + 1;

        U_cur = U_pred + delta_U;
        V_cur = V_pred + gamma_nm*dt*A_new;

        % Superstructure restoring force
        F_restore = K_global * U_cur;

        % Damping force
        F_damp = C_damp * V_cur;

        % Inertial force
        F_inertia = M_global * A_new;

        % External force from ground motion (effective earthquake force)
        F_ext = -M_global * (ag * ones(n_dofs, 1));

        % Isolator forces and tangent stiffness
        F_iso_vec = zeros(n_dofs, 1);
        K_iso_tan = zeros(n_dofs, n_dofs);

        for i_iso = 1:n_isolators
            node = i_iso;  % Level 0 node
            dof_ux = 3*(node-1) + 1;
            dof_uy = 3*(node-1) + 2;

            u_iso = U_cur(dof_ux);
            du_iso = U_cur(dof_ux) - U_prev(dof_ux);

            % Bouc-Wen state update: iterate z to convergence
            z_old = z_bw(i_iso);
            z_iter = z_old;
            for bw_iter = 1:30
                abs_z = abs(z_iter);
                % Residual: z_new - z_old - Ao*du + beta*|du|*|z|*z + gamma*du*z^2
                R_bw = z_iter - z_old - Ao*du_iso + beta*abs(du_iso)*abs_z*z_iter + gamma_bw*du_iso*z_iter^2;
                % Jacobian: dR/dz = 1 + 2*beta*|du|*|z| + 2*gamma*du*z
                J_bw = 1 + 2*beta*abs(du_iso)*abs_z + 2*gamma_bw*du_iso*z_iter;
                dz_corr = -R_bw / J_bw;
                z_iter = z_iter + dz_corr;
                if abs(dz_corr) < 1e-14
                    break;
                end
            end
            z_new = z_iter;

            % Isolator horizontal force
            F_h = alpha * k0 * u_iso + (1 - alpha) * k0 * z_new;

            % Tangent stiffness: analytic derivative
            % K_tan = alpha*k0 + (1-alpha)*k0 * dz_new/du
            % dz/du = (Ao - gamma*|z|^n) / (1 + 2*beta*|du|*|z| + 2*gamma*du*z)
            % But for Newton, use the elastic stiffness for robustness
            k_tan_h = k0;

            F_iso_vec(dof_ux) = F_iso_vec(dof_ux) + F_h;
            K_iso_tan(dof_ux, dof_ux) = K_iso_tan(dof_ux, dof_ux) + k_tan_h;

            % Vertical stiffness
            F_iso_vec(dof_uy) = F_iso_vec(dof_uy) + kv * U_cur(dof_uy);
            K_iso_tan(dof_uy, dof_uy) = K_iso_tan(dof_uy, dof_uy) + kv;
        end

        % Residual
        Residual = F_ext - F_inertia - F_damp - F_restore - F_iso_vec;

        % Effective stiffness with isolator contribution
        K_eff = K_eff_linear + K_iso_tan;

        % Solve for displacement correction
        dU_solve = K_eff \ Residual;

        delta_U = delta_U + dU_solve;
        A_new = A_new + a0 * dU_solve;

        % Convergence check (relative criterion)
        F_restore_total = abs(F_restore + F_iso_vec);
        R_norm = norm(Residual, inf);
        denom = max([1.0, norm(F_ext, inf), norm(F_restore + F_iso_vec, inf)]);
        rel_err = R_norm / denom;

        if rel_err < tol_nr
            converged = true;
        end
    end

    if converged
        % Accept the solution
        U = U_pred + delta_U;
        V = V_pred + gamma_nm * dt * A_new;
        A_acc = A_new;

        % Update Bouc-Wen state at converged displacement (use total incremental displ)
        for i_iso = 1:n_isolators
            node = i_iso;
            dof_ux = 3*(node-1) + 1;
            du_iso_final = U(dof_ux) - U_prev(dof_ux);

            z_old = z_bw(i_iso);
            z_iter = z_old;
            for bw_iter = 1:30
                abs_z = abs(z_iter);
                R_bw = z_iter - z_old - Ao*du_iso_final + ...
                    beta*abs(du_iso_final)*abs_z*z_iter + gamma_bw*du_iso_final*z_iter^2;
                J_bw = 1 + 2*beta*abs(du_iso_final)*abs_z + 2*gamma_bw*du_iso_final*z_iter;
                dz_corr = -R_bw / J_bw;
                z_iter = z_iter + dz_corr;
                if abs(dz_corr) < 1e-14
                    break;
                end
            end
            z_bw(i_iso) = z_iter;
        end
    else
        failed_steps = failed_steps + 1;
    end

    % Store results
    store_t(step) = step * dt;

    % Roof displacement: level 3, center node (node 11)
    roof_node = 3*(n_bays+1) + 2;  % = 11
    roof_dof_ux = 3*(roof_node-1) + 1;  % = 31
    store_roof_ux(step) = U(roof_dof_ux);

    % Isolator data and base shear
    base_shear_total = 0;
    for i_iso = 1:n_isolators
        node = i_iso;
        dof_ux = 3*(node-1) + 1;
        u_iso = U(dof_ux);
        F_iso = alpha * k0 * u_iso + (1-alpha) * k0 * z_bw(i_iso);

        store_iso_ux(step, i_iso) = u_iso;
        store_iso_F(step, i_iso) = F_iso;
        store_iso_z(step, i_iso) = z_bw(i_iso);
        base_shear_total = base_shear_total + F_iso;
    end
    store_base_shear(step) = base_shear_total;

    % Story drifts (from center column line)
    ux_center = zeros(4,1);
    for lv = 0:3
        node_center = lv*(n_bays+1) + 2;  % Center column line (column 2)
        dof = 3*(node_center-1) + 1;
        ux_center(lv+1) = U(dof);
    end
    drift1 = (ux_center(2) - ux_center(1)) / H_story;
    drift2 = (ux_center(3) - ux_center(2)) / H_story;
    drift3 = (ux_center(4) - ux_center(3)) / H_story;

    store_drift1(step) = drift1;
    store_drift2(step) = drift2;
    store_drift3(step) = drift3;

    % Update peaks
    peak_roof_ux = max(peak_roof_ux, abs(store_roof_ux(step)));
    peak_iso_ux  = max(peak_iso_ux, max(abs(store_iso_ux(step,:))));
    peak_base_shear = max(peak_base_shear, abs(base_shear_total));
    peak_drift1 = max(peak_drift1, abs(drift1));
    peak_drift2 = max(peak_drift2, abs(drift2));
    peak_drift3 = max(peak_drift3, abs(drift3));

    store_n_iter(step) = n_iter;
    store_converged(step) = converged;

    if mod(step, 200) == 0
        fprintf('Step %4d/%d: t=%.1fs, roof=%.4fm, shear=%.1fkN, iter=%d\n', ...
            step, n_steps, step*dt, store_roof_ux(step), base_shear_total/1000, n_iter);
    end
end

fprintf('Analysis complete. Failed: %d/%d steps\n', failed_steps, n_steps);

%% ========================================================================
% 16. OUTPUT SUMMARY
%  ========================================================================
fprintf('\n=== Peak Response Summary ===\n');
fprintf('Fundamental period (elastic): %.4f s\n', T_periods(1));
fprintf('Peak roof displacement:       %.4f m\n', peak_roof_ux);
fprintf('Peak isolation displacement:  %.4f m\n', peak_iso_ux);
fprintf('Peak total base shear:        %.2f kN\n', peak_base_shear/1000);
fprintf('Max drift ratio story 1:      %.6f\n', peak_drift1);
fprintf('Max drift ratio story 2:      %.6f\n', peak_drift2);
fprintf('Max drift ratio story 3:      %.6f\n', peak_drift3);
fprintf('Failed steps:                 %d\n', failed_steps);

%% ========================================================================
% 17. SAVE OUTPUTS
%  ========================================================================
response = struct();
response.software = 'MATLAB';
response.fundamental_period_s = T_periods(1);
response.periods_s = T_periods';
response.frequencies_hz = freqs';
response.rayleigh_alpha_m = alpha_m;
response.rayleigh_beta_k = beta_k;
response.peak_roof_displacement_m = peak_roof_ux;
response.peak_isolation_displacement_m = peak_iso_ux;
response.peak_base_shear_N = peak_base_shear;
response.peak_drift_ratio_story1 = peak_drift1;
response.peak_drift_ratio_story2 = peak_drift2;
response.peak_drift_ratio_story3 = peak_drift3;
response.failed_steps = failed_steps;
response.total_steps = n_steps;
response.convergence_tolerance = tol_nr;
response.integration_scheme = 'Newmark-beta_average';
response.bouc_wen_update = 'backward_Euler_with_Newton_subiteration';
response.bouc_wen_equation = 'dz = Ao*du - (beta*|du|*|z|^(n-1)*z + gamma*du*|z|^n)';
response.n_steps = n_steps;
response.dt_s = dt;
response.pga_target_m_s2 = target_pga;
response.model_type = '2D_elastic_frame_with_BoucWen_isolators';
response.diaphragm_method = 'penalty_K1e14';

fid = fopen('response.json', 'w');
fprintf(fid, '%s', jsonencode(response, 'PrettyPrint', true));
fclose(fid);

% Time history file
th_header = 'time roof_ux iso_ux1 iso_ux2 iso_ux3 iso_F1 iso_F2 iso_F3 base_shear drift1 drift2 drift3 iso_z1 iso_z2 iso_z3 n_iter';
th_data = [store_t, store_roof_ux, store_iso_ux, store_iso_F, ...
    store_base_shear, store_drift1, store_drift2, store_drift3, ...
    store_iso_z, store_n_iter];

fid = fopen('response_time_history.txt', 'w');
fprintf(fid, '# %s\n', th_header);
fclose(fid);
dlmwrite('response_time_history.txt', th_data, '-append', 'delimiter', '\t', 'precision', '%.12e');

% Individual isolator hysteresis files
for i_iso = 1:n_isolators
    iso_data = [store_t, store_iso_ux(:,i_iso), store_iso_F(:,i_iso), store_iso_z(:,i_iso)];
    fname = sprintf('isolator_%d_hysteresis.txt', i_iso);
    fid = fopen(fname, 'w');
    fprintf(fid, '# time ux F_h z\n');
    fclose(fid);
    dlmwrite(fname, iso_data, '-append', 'delimiter', '\t', 'precision', '%.12e');
end

fprintf('\n=== Output files written ===\n');
fprintf('  response.json\n  response_time_history.txt\n  isolator_*_hysteresis.txt\n');

%% ========================================================================
% 18. PLOTS
%  ========================================================================
figure(1); clf;
subplot(3,2,1);
plot(store_t, store_roof_ux); xlabel('Time (s)'); ylabel('Roof UX (m)');
title('Roof Displacement'); grid on;

subplot(3,2,2);
plot(store_t, store_base_shear/1000); xlabel('Time (s)'); ylabel('Base Shear (kN)');
title('Total Base Shear'); grid on;

subplot(3,2,3);
plot(store_t, store_iso_ux); xlabel('Time (s)'); ylabel('Iso UX (m)');
title('Isolation Displacement'); grid on; legend('Left','Center','Right');

subplot(3,2,4);
plot(store_iso_ux(:,2), store_iso_F(:,2)/1000);
xlabel('Iso UX (m)'); ylabel('Iso Force (kN)');
title('Center Isolator Hysteresis'); grid on;

subplot(3,2,5);
plot(store_t, [store_drift1, store_drift2, store_drift3]);
xlabel('Time (s)'); ylabel('Drift Ratio');
title('Story Drift Ratios'); grid on; legend('S1','S2','S3');

subplot(3,2,6);
plot(store_t, store_iso_z); xlabel('Time (s)'); ylabel('BW state z');
title('Bouc-Wen State'); grid on; legend('Left','Center','Right');

saveas(gcf, 'response_plots.png');
fprintf('  response_plots.png\n');

fprintf('\n=== MATLAB analysis complete ===\n');
