% Seismic isolation 4-story building analysis
% Newmark-beta average acceleration method
clear; close all;

%% Load ground motion
gm = load('../../input_data/GM1.txt');
dt = 0.02;          % time step, s
g = 9.81;           % gravity, m/s^2
target_pga = 0.40 * g;
max_abs_gm = max(abs(gm));
scale = target_pga / max_abs_gm;
ag = gm * scale;    % normalized ground acceleration, m/s^2
N = length(ag);
t = (0:N-1)' * dt;
fprintf('GM1: N=%d, max_abs=%.6f, scale=%.6f, target_PGA=%.4f\n', N, max_abs_gm, scale, target_pga);

%% System matrices (4 DOF: 1=isolation, 2=story2, 3=story3, 4=story4/top)
M = diag([500, 1000, 1000, 1000]);  % kg

k_iso = 200e3;   % N/m
k_story = 500e3; % N/m

K = [k_iso+k_story, -k_story,        0,         0;
     -k_story,     2*k_story, -k_story,         0;
      0,          -k_story,  2*k_story, -k_story;
      0,             0,     -k_story,   k_story];

%% Natural frequencies and Rayleigh damping
[phi, omega2] = eig(K, M);
omega = sqrt(diag(omega2));
omega = sort(omega);
omega1 = omega(1);
omega2 = omega(2);
fprintf('Natural frequencies (rad/s): %.4f, %.4f, %.4f, %.4f\n', omega);
fprintf('Natural frequencies (Hz): %.4f, %.4f, %.4f, %.4f\n', omega/(2*pi));

xi = 0.05;
alpha = 2*xi*omega1*omega2/(omega1+omega2);
beta_ray = 2*xi/(omega1+omega2);
fprintf('Rayleigh: alpha=%.6f, beta=%.6e\n', alpha, beta_ray);

% Rayleigh damping matrix
C_R = alpha * M + beta_ray * K;

% Add isolation-layer damping coefficient to C(1,1)
c_iso = 2000;  % N*s/m
C = C_R;
C(1,1) = C(1,1) + c_iso;
fprintf('Added isolation damping c=%g to C(1,1)\n', c_iso);

%% Newmark-beta time integration (average acceleration)
gamma_nm = 0.5;
beta_nm = 0.25;

ndof = 4;
u = zeros(ndof, N);    % relative displacement
v = zeros(ndof, N);    % relative velocity
a = zeros(ndof, N);    % relative acceleration

% Initial acceleration (from equilibrium at t=0)
influence = ones(ndof, 1);
a(:,1) = M \ (-C*v(:,1) - K*u(:,1) - M*influence*ag(1));

% Effective stiffness (constant for linear system)
K_eff = M + gamma_nm*dt*C + beta_nm*dt^2*K;

fprintf('Newmark integration: %d steps, dt=%.3f s\n', N-1, dt);

% Time-stepping loop
for n = 1:N-1
    % Predictor step
    u_pred = u(:,n) + dt*v(:,n) + 0.5*dt^2*(1-2*beta_nm)*a(:,n);
    v_pred = v(:,n) + dt*(1-gamma_nm)*a(:,n);

    % Effective force
    F_eff = -M*influence*ag(n+1) - C*v_pred - K*u_pred;

    % Solve for new acceleration
    a(:,n+1) = K_eff \ F_eff;

    % Corrector step
    u(:,n+1) = u_pred + beta_nm*dt^2*a(:,n+1);
    v(:,n+1) = v_pred + gamma_nm*dt*a(:,n+1);
end

%% Extract results
top_dis = u(4,:)';          % top-story relative displacement
top_acc_rel = a(4,:)';          % top-story relative acceleration
top_acc_abs = top_acc_rel + ag;  % top-story absolute acceleration
iso_dis = u(1,:)';               % isolation-layer relative displacement

fprintf('Peak top displacement: %.6f m\n', max(abs(top_dis)));
fprintf('Peak top rel acceleration: %.6f m/s^2\n', max(abs(top_acc_rel)));
fprintf('Peak top abs acceleration: %.6f m/s^2 (%.4f g)\n', max(abs(top_acc_abs)), max(abs(top_acc_abs))/g);
fprintf('Peak isolation displacement: %.6f m\n', max(abs(iso_dis)));

%% Write output files
% topStoDis.txt: time + top displacement
fid = fopen('topStoDis.txt', 'w');
fprintf(fid, '%.10e %.10e\n', [t'; top_dis']);
fclose(fid);

% topStoAcc.txt: time + top absolute acceleration
fid = fopen('topStoAcc.txt', 'w');
fprintf(fid, '%.10e %.10e\n', [t'; top_acc_rel']);
fclose(fid);

% isolation displacement
fid = fopen('isoDis.txt', 'w');
fprintf(fid, '%.10e %.10e\n', [t'; iso_dis']);
fclose(fid);

% Response JSON
resp = struct();
resp.time = t;
resp.top_displacement = top_dis;
resp.top_rel_acceleration = top_acc_rel;
resp.top_abs_acceleration = top_acc_abs;
resp.isolation_displacement = iso_dis;
resp.input_pga = target_pga;
resp.dt = dt;
resp.n_samples = N;
resp.peak_top_dis = max(abs(top_dis));
resp.peak_top_rel_acc = max(abs(top_acc_rel));
resp.peak_top_abs_acc = max(abs(top_acc_abs));
resp.peak_iso_dis = max(abs(iso_dis));
resp.omega = omega;
resp.alpha = alpha;
resp.beta = beta_ray;
resp.c_iso_added = c_iso;

fid = fopen('response.json', 'w');
fprintf(fid, '%s', jsonencode(resp, 'PrettyPrint', true));
fclose(fid);

fprintf('MATLAB analysis complete. Output files written.\n');
