% Nonlinear MDOF time-history analysis for 4-story isolated structure
% Upper 3 stories: ideal elastic-plastic (b=0)
% Isolation layer: linear (keq + damping c)
% Rayleigh damping ratio = 0.05
% Newmark-beta with Newton-Raphson iteration

clear; close all;

%% Parameters from supplementary/Example4_2.m
mass_isolation = 250; % t
mass_stories = [270, 270, 180]; % t
k_stories = [245, 195, 98] * 1e6; % N/m  (k2, k3, k4)
fy_stories = [1225, 975, 490] * 1e3; % N  (fy2, fy3, fy4)
keq = 50 * 1e6; % N/m  (isolation stiffness)
c_iso = 1000 * 1e3; % N/(m/s)  (isolation damping coefficient)
damping_ratio = 0.05;
b_ratio = [0, 0, 0]; % ideal elastic-plastic
dt = 0.01; % s
g = 9.81; % m/s^2
PGA = 0.40 * g;

%% Load and normalize ground motion
accel_data = load('../../input_data/Northridge_01_NO_968.txt');
accel_data = accel_data / max(abs(accel_data)) * PGA;
n_steps = length(accel_data);
time = (0:n_steps-1)' * dt;

%% System matrices
n_dof = 4;
n_story = 3; % upper 3 stories (DOFs 2,3,4)
m_all = [mass_isolation, mass_stories] * 1e3; % kg
M = diag(m_all);

% Initial elastic stiffness matrix
K0 = zeros(n_dof);
K0(1,1) = keq + k_stories(1);
K0(1,2) = -k_stories(1);
K0(2,1) = -k_stories(1);
for i = 2:n_dof-1
    K0(i,i) = k_stories(i-1) + k_stories(i);
    K0(i,i-1) = -k_stories(i-1);
    K0(i,i+1) = -k_stories(i);
end
K0(n_dof,n_dof) = k_stories(end);
K0(n_dof,n_dof-1) = -k_stories(end);

% Rayleigh damping
[V, D] = eig(K0, M);
omega = sqrt(diag(D));
omega1 = omega(1);
omega2 = omega(2);
alpha = 2*omega1*omega2*(damping_ratio*omega2 - damping_ratio*omega1)/(omega2^2 - omega1^2);
beta_r = 2*(damping_ratio*omega2 - damping_ratio*omega1)/(omega2^2 - omega1^2);
C = alpha*M + beta_r*K0;
C(1,1) = C(1,1) + c_iso;

%% Newmark-beta constants
beta_n = 1/4;  % average acceleration
gamma_n = 1/2;
c0 = 1/(beta_n*dt^2);
c1 = gamma_n/(beta_n*dt);
c2 = 1/(beta_n*dt);
c3 = 1/(2*beta_n) - 1;
c4 = gamma_n/beta_n - 1;
c5 = dt*(gamma_n/(2*beta_n) - 1);

%% Constant part of effective stiffness: K_const = c0*M + c1*C
K_const = c0*M + c1*C;

%% Initial conditions
u_prev = zeros(n_dof, 1);
v_prev = zeros(n_dof, 1);
a_prev = -ones(n_dof, 1) * accel_data(1);

%% Nonlinear state variables (for stories 2, 3, 4)
lastShear = zeros(n_story, 1);
lastRelDisp = zeros(n_story, 1);

%% Storage arrays
disp_top = zeros(n_steps, 1);
acc_top = zeros(n_steps, 1);
disp_iso = zeros(n_steps, 1);
base_shear = zeros(n_steps, 1);
yield_flag = zeros(n_steps, n_story);
convergence_iter = zeros(n_steps, 1);
failed_steps = 0;

%% Newton-Raphson parameters
max_iter = 50;
tol_res = 1e-6;

%% Time-stepping loop
for step = 1:n_steps
    ag = accel_data(step);

    % Effective force predictor (from previous converged state)
    F_eff = -M * ones(n_dof,1) * ag ...
            + M * (c0*u_prev + c2*v_prev + c3*a_prev) ...
            + C * (c1*u_prev + c4*v_prev + c5*a_prev);

    % Initial guess: previous converged displacement
    u_cur = u_prev;
    relDisp_cur = zeros(n_story, 1);
    relDisp_cur(1) = u_cur(2) - u_cur(1);
    relDisp_cur(2) = u_cur(3) - u_cur(2);
    relDisp_cur(3) = u_cur(4) - u_cur(3);

    converged = false;
    for iter = 1:max_iter
        % Compute nonlinear story forces and tangent stiffness
        kt_story = zeros(n_story, 1);
        shear = zeros(n_story, 1);

        for s = 1:n_story
            deltaDisp = relDisp_cur(s) - lastRelDisp(s);
            trialShear = lastShear(s) + k_stories(s) * deltaDisp;

            % Ideal elastic-plastic (b=0): yield bounds [-fy, +fy]
            fy_s = fy_stories(s);
            shear(s) = max(-fy_s, min(fy_s, trialShear));

            if abs(shear(s) - trialShear) < 1e-3
                kt_story(s) = k_stories(s);
            else
                kt_story(s) = 0;
            end
        end

        % Track yielding
        for s = 1:n_story
            if abs(shear(s)) >= fy_stories(s) * 0.999
                yield_flag(step, s) = 1;
            end
        end

        % Build restoring force vector R(u)
        R = zeros(n_dof, 1);
        R(1) = keq * u_cur(1) - shear(1);
        R(2) = shear(1) - shear(2);
        R(3) = shear(2) - shear(3);
        R(4) = shear(3);

        % Build tangent stiffness matrix
        Kt = zeros(n_dof);
        Kt(1,1) = keq + kt_story(1);
        Kt(1,2) = -kt_story(1);
        for s = 1:n_story
            if s == 1
                Kt(2,1) = -kt_story(1);
                Kt(2,2) = kt_story(1) + kt_story(2);
            elseif s == 2
                Kt(2,3) = -kt_story(2);
                Kt(3,2) = -kt_story(2);
                Kt(3,3) = kt_story(2) + kt_story(3);
            elseif s == 3
                Kt(3,4) = -kt_story(3);
                Kt(4,3) = -kt_story(3);
                Kt(4,4) = kt_story(3);
            end
        end

        % Residual: g(u) = K_const*u + R(u) - F_eff
        res = K_const * u_cur + R - F_eff;

        if norm(res, inf) < tol_res
            converged = true;
            break;
        end

        % Jacobian: J = K_const + Kt
        J = K_const + Kt;
        du = J \ (-res);
        u_cur = u_cur + du;

        % Update interstory drifts for new u_cur
        relDisp_cur(1) = u_cur(2) - u_cur(1);
        relDisp_cur(2) = u_cur(3) - u_cur(2);
        relDisp_cur(3) = u_cur(4) - u_cur(3);
    end

    if ~converged
        failed_steps = failed_steps + 1;
    end
    convergence_iter(step) = iter;

    % Compute new acceleration and velocity
    a_cur = c0*(u_cur - u_prev) - c2*v_prev - c3*a_prev;
    v_cur = v_prev + dt*((1-gamma_n)*a_prev + gamma_n*a_cur);

    % Store results
    disp_top(step) = u_cur(4);
    acc_top(step) = a_cur(4);
    disp_iso(step) = u_cur(1);
    base_shear(step) = keq * u_cur(1);

    % Update nonlinear state for next step
    relDisp_final = zeros(n_story, 1);
    relDisp_final(1) = u_cur(2) - u_cur(1);
    relDisp_final(2) = u_cur(3) - u_cur(2);
    relDisp_final(3) = u_cur(4) - u_cur(3);

    % Recompute shear at converged state
    for s = 1:n_story
        deltaDisp = relDisp_final(s) - lastRelDisp(s);
        trialShear = lastShear(s) + k_stories(s) * deltaDisp;
        fy_s = fy_stories(s);
        lastShear(s) = max(-fy_s, min(fy_s, trialShear));
    end
    lastRelDisp = relDisp_final;

    % Advance state
    u_prev = u_cur;
    v_prev = v_cur;
    a_prev = a_cur;
end

%% Save results
topStoDisIso1 = [time, disp_top];
topStoAccIso1 = [time, acc_top];
save('topStoDisIso1.txt', 'topStoDisIso1', '-ascii', '-double');
save('topStoAccIso1.txt', 'topStoAccIso1', '-ascii', '-double');

% Save isolation displacement
iso_disp_out = [time, disp_iso];
save('iso_displacement.txt', 'iso_disp_out', '-ascii', '-double');

% Save hysteresis data for story 2
hyst_data = [time, base_shear, disp_iso];
save('hysteresis_layer1.txt', 'hyst_data', '-ascii', '-double');

%% Print summary
fprintf('Analysis complete.\n');
fprintf('Steps: %d, Failed: %d\n', n_steps, failed_steps);
fprintf('Max top displacement: %.6f m\n', max(abs(disp_top)));
fprintf('Max top acceleration: %.6f m/s^2\n', max(abs(acc_top)));
fprintf('Max isolation displacement: %.6f m\n', max(abs(disp_iso)));
fprintf('Yielding steps - Story2: %d, Story3: %d, Story4: %d\n', ...
    sum(yield_flag(:,1)), sum(yield_flag(:,2)), sum(yield_flag(:,3)));

%% Plot
figure;
subplot(3,1,1);
plot(time, disp_top);
title('Top Story Displacement'); xlabel('Time (s)'); ylabel('Disp (m)');

subplot(3,1,2);
plot(time, acc_top);
title('Top Story Acceleration'); xlabel('Time (s)'); ylabel('Acc (m/s^2)');

subplot(3,1,3);
plot(time, disp_iso);
title('Isolation Displacement'); xlabel('Time (s)'); ylabel('Disp (m)');

saveas(gcf, 'response_plot.png');
