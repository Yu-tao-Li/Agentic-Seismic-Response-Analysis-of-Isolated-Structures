% Load seismic data
data = load('Northridge_01_NO_968.txt');
dt = 0.01; % Time interval in seconds
g = 9.81; % Acceleration due to gravity in m/s^2
PGA = 0.40 * g; % Peak Ground Acceleration in m/s^2

% Normalize the seismic wave data
acceleration_data = data(:,1);
time = (0:length(acceleration_data)-1)' * dt;
acceleration_data = acceleration_data / max(abs(acceleration_data)) * PGA;

% Structure parameters
m = [250 270 270 180] * 1e3; % Mass of each story in kg
k = [50 245 195 98] * 1e6; % Stiffness of each story in N/m
fy = [500 1225 975 490] * 1e3; % Yield capacity of each story in N
b = [0.05 0 0 0]; % Post-yield stiffness hardening coefficient

% Rayleigh damping parameters
damping_ratio = 0.05;

% Mass matrix
M = diag(m);

% Stiffness matrix
K = diag([k(1)+k(2), k(2)+k(3), k(3)+k(4), k(4)]) - diag(k(2:end), 1) - diag(k(2:end), -1);

% Compute eigenvalues and eigenvectors
[V, D] = eig(K, M);
omega = sqrt(diag(D));

% Rayleigh damping coefficients
alpha = 2 * damping_ratio * omega(1) * omega(2) / (omega(1) + omega(2));
beta = 2 * damping_ratio / (omega(1) + omega(2));

% Damping matrix
C = alpha * M + beta * K;

% Newmark-beta method parameters
beta_n = 1/4;
gamma = 1/2;

% Initial conditions
u = zeros(4, 1);
u_dot = zeros(4, 1);
u_dot_dot = zeros(4, 1);

% Effective stiffness matrix
K_eff = K + gamma / (beta_n * dt) * C + 1 / (beta_n * dt^2) * M;

% Time-stepping loop
displacement = zeros(length(time), 1);
acceleration = zeros(length(time), 1);

for i = 1:length(time)
    % Calculate effective force vector
    F_eff = -M * acceleration_data(i) * ones(4, 1) + M * (1 / (beta_n * dt^2) * u + 1 / (beta_n * dt) * u_dot + (1 / (2 * beta_n) - 1) * u_dot_dot) + C * (gamma / (beta_n * dt) * u + (gamma / beta_n - 1) * u_dot + dt * (gamma / (2 * beta_n) - 1) * u_dot_dot);
    
    % Solve for displacement
    u_new = K_eff \ F_eff;
    
    % Calculate new acceleration and velocity
    u_dot_dot_new = 1 / (beta_n * dt^2) * (u_new - u) - 1 / (beta_n * dt) * u_dot - (1 / (2 * beta_n) - 1) * u_dot_dot;
    u_dot_new = u_dot + dt * ((1 - gamma) * u_dot_dot + gamma * u_dot_dot_new);
    
    % Update displacements, velocities, and accelerations
    u = u_new;
    u_dot = u_dot_new;
    u_dot_dot = u_dot_dot_new;
    
    % Store top story displacement and acceleration
    displacement(i) = u(end);
    acceleration(i) = u_dot_dot(end);
end

% Save results to files
topStoDisIso2 = [time, displacement];
topStoAccIso2 = [time, acceleration];
save('topStoDisIso2.txt', 'topStoDisIso2', '-ascii');
save('topStoAccIso2.txt', 'topStoAccIso2', '-ascii');

% Plot results
figure;
subplot(2, 1, 1);
plot(time, displacement);
title('Top Story Displacement Time History');
xlabel('Time (s)');
ylabel('Displacement (m)');

subplot(2, 1, 2);
plot(time, acceleration);
title('Top Story Acceleration Time History');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
