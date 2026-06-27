% Linear MDOF analysis for direct comparison with OpenSeesPy
clear;

%% Parameters
mass_isolation = 250; % t
mass_stories = [270, 270, 180]; % t
k_stories = [245, 195, 98] * 1e6;
keq = 50 * 1e6;
c_iso = 1000 * 1e3;
damping_ratio = 0.05;
dt = 0.01;
g = 9.81;
PGA = 0.40 * g;

accel_data = load('../../input_data/Northridge_01_NO_968.txt');
accel_data = accel_data / max(abs(accel_data)) * PGA;
n_steps = length(accel_data);
time = (0:n_steps-1)' * dt;

n_dof = 4;
m_all = [mass_isolation, mass_stories] * 1e3;
M = diag(m_all);

K0 = zeros(n_dof);
K0(1,1) = keq + k_stories(1);
K0(1,2) = -k_stories(1); K0(2,1) = -k_stories(1);
for i = 2:n_dof-1
    K0(i,i) = k_stories(i-1) + k_stories(i);
    K0(i,i-1) = -k_stories(i-1); K0(i,i+1) = -k_stories(i);
end
K0(n_dof,n_dof) = k_stories(end);
K0(n_dof,n_dof-1) = -k_stories(end);

[V, D] = eig(K0, M);
omega = sqrt(diag(D));
omega1 = omega(1); omega2 = omega(2);
alpha = 2*omega1*omega2*(damping_ratio*omega2 - damping_ratio*omega1)/(omega2^2 - omega1^2);
beta = 2*(damping_ratio*omega2 - damping_ratio*omega1)/(omega2^2 - omega1^2);

C = alpha*M + beta*K0;
C(1,1) = C(1,1) + c_iso;

beta_n = 1/4; gamma_n = 1/2;
c0 = 1/(beta_n*dt^2);
c1 = gamma_n/(beta_n*dt);
c2 = 1/(beta_n*dt);
c3 = 1/(2*beta_n) - 1;
c4 = gamma_n/beta_n - 1;
c5 = dt*(gamma_n/(2*beta_n) - 1);

K_eff = K0 + c1*C + c0*M;

u = zeros(n_dof,1); v = zeros(n_dof,1); a = -ones(n_dof,1)*accel_data(1);

dmax = 0; amax = 0;
for i = 1:n_steps
    ag = accel_data(i);
    F_eff = -M*ones(n_dof,1)*ag ...
            + M*(c0*u + c2*v + c3*a) ...
            + C*(c1*u + c4*v + c5*a);
    u_new = K_eff \ F_eff;
    a_new = c0*(u_new - u) - c2*v - c3*a;
    v_new = v + dt*((1-gamma_n)*a + gamma_n*a_new);

    u = u_new; v = v_new; a = a_new;
    dmax = max(dmax, abs(u(4)));
    amax = max(amax, abs(a(4)));
end

fprintf('MATLAB linear: d_max=%.6f, a_max=%.6f\n', dmax, amax);
