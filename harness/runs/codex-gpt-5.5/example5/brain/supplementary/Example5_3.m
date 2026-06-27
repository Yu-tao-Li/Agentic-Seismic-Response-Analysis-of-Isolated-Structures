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
beta_n = 0.25;
gamma = 0.5;

% Initial conditions
u = zeros(4, 1);
u_dot = zeros(4, 1);
u_dot_dot = M \ (-acceleration_data(1) * diag(M) - C * u_dot);

% Initialize arrays to store results
displacement = zeros(length(time), 1);
acceleration = zeros(length(time), 1);
displacement(1) = u(end);
acceleration(1) = u_dot_dot(end);

% Initial restoring force
Rs = zeros(4, 1);

for j = 1:length(time)-1
    % External load increment
    delta_P = -(acceleration_data(j+1) - acceleration_data(j)) * diag(M);
    
    % Effective load increment
    delta_P_eff = delta_P + ((1 / (beta_n * dt)) * M + (gamma / beta_n) * C) * u_dot + ((1 / (2 * beta_n)) * M + dt * (gamma / (2 * beta_n) - 1) * C) * u_dot_dot;
    
    % Effective stiffness matrix
    K_eff = K + (1 / (beta_n * dt^2)) * M + (gamma / (beta_n * dt)) * C;

    % Newton-Raphson iteration
    u_trial = u;
    Rs_trial = Rs;
    delta_R_trial = delta_P_eff;
    
    for m = 1:100 % maximum number of iterations
        if sum(abs(delta_R_trial)) < 1e-5
            break;
        end
        
        delta_u = K_eff \ delta_R_trial;
        u_trial_n = u_trial + delta_u;
        
        % Update restoring force and stiffness
        [Rs_trial_n, k_trial] = computeRsStiffness(Rs_trial, u_trial, u_trial_n, fy, k , b);
        K_trial_eff = computeStiffnessMatrix(k_trial) + (1 / (beta_n * dt^2)) * M + (gamma / (beta_n * dt)) * C;
        
        delta_F = Rs_trial_n - Rs_trial + (gamma / (beta_n * dt)) * C * delta_u + (1 / (beta_n * dt^2)) * M * delta_u;
        delta_R_trial = delta_R_trial - delta_F;
        
        % Update trial values
        u_trial = u_trial_n;
        Rs_trial = Rs_trial_n;
        K_eff = K_trial_eff;
    end
    
    % Update displacement, velocity, and acceleration
    delta_u = u_trial - u;
    delta_u_dot = (gamma / (beta_n * dt)) * delta_u - (gamma / beta_n) * u_dot + (1 - gamma / (2 * beta_n)) * dt * u_dot_dot;
    delta_u_dot_dot = (1 / (beta_n * dt^2)) * delta_u - (1 / (beta_n * dt)) * u_dot - (1 / (2 * beta_n)) * u_dot_dot;
    
    u = u_trial;
    u_dot = u_dot + delta_u_dot;
    u_dot_dot = u_dot_dot + delta_u_dot_dot;
    
    % Store top story displacement and acceleration
    displacement(j+1) = u(end);
    acceleration(j+1) = u_dot_dot(end);
    
    % Update restoring force
    Rs = Rs_trial;
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
