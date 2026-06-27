% Define the parameters for the isolation layer and stories
mass_isolation = 500; % kg
stiffness_isolation = 200 * 1e3; % N/m
damping_isolation = 2; % kN/(m·s)
damping_isolation = damping_isolation * 1e3; % Convert to N/(m·s)

mass_story = 1000; % kg
stiffness_story = 500 * 1e3; % N/m
damping_ratio = 0.05; % Rayleigh damping ratio
dt = 0.02; % time interval in seconds
g = 9.81; % gravity in m/s^2
PGA = 0.40 * g; % Peak Ground Acceleration in m/s^2

% Read the seismic wave data
acceleration_data = load('GM1.txt');
time = (0:length(acceleration_data)-1)' * dt;

% Normalize the acceleration data to the specified PGA
acceleration_data = acceleration_data / max(abs(acceleration_data)) * PGA;

% Define the number of stories including the isolation layer
n_stories = 4;

% Mass matrix
M = mass_story * eye(n_stories);
M(1,1) = mass_isolation; % Mass of the isolation layer

% Stiffness matrix
K = zeros(n_stories);
K(1,1) = stiffness_isolation + stiffness_story;
K(1,2) = -stiffness_story;
K(2,1) = -stiffness_story;
for i = 2:n_stories-1
    K(i,i) = 2 * stiffness_story;
    K(i,i-1) = -stiffness_story;
    K(i,i+1) = -stiffness_story;
end
K(n_stories,n_stories) = stiffness_story;
K(n_stories,n_stories-1) = -stiffness_story;

% Compute the eigenvalues and eigenvectors
[V, D] = eig(K, M);
omega = sqrt(diag(D));

% Rayleigh damping coefficients
omega1 = omega(1);
omega2 = omega(2);
alpha = 2 * omega1 * omega2 * (damping_ratio * omega2 - damping_ratio * omega1) / (omega2^2 - omega1^2);
beta = 2 * (damping_ratio * omega2 - damping_ratio * omega1) / (omega2^2 - omega1^2);

% Damping matrix
C = alpha * M + beta * K;
C(1,1) = C(1,1) + damping_isolation; % Adding damping coefficient of the isolation layer

% Newmark-beta method parameters
beta_n = 1/4;
gamma = 1/2;

% Initial conditions
u = zeros(n_stories, 1);
u_dot = zeros(n_stories, 1);
u_dot_dot = zeros(n_stories, 1);
u_dot_dot_p = zeros(n_stories, 1);

% Initialize arrays to store results
displacement = zeros(length(time), 1);
acceleration = zeros(length(time), 1);

% Effective stiffness matrix
K_eff = K + gamma / (beta_n * dt) * C + 1 / (beta_n * dt^2) * M;

% Time-stepping loop
for i = 1:length(time)
    % Calculate effective force vector
    F_eff = -M * acceleration_data(i) * ones(n_stories, 1) + M * (1 / (beta_n * dt^2) * u + 1 / (beta_n * dt) * u_dot + (1 / (2 * beta_n) - 1) * u_dot_dot) + C * (gamma / (beta_n * dt) * u + (gamma / beta_n - 1) * u_dot + dt * (gamma / (2 * beta_n) - 1) * u_dot_dot);
    
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
topFloDisIso = [time, displacement];
topFloAccIso = [time, acceleration];
save('topFloDisIso.txt', 'topFloDisIso', '-ascii');
save('topFloAccIso.txt', 'topFloAccIso', '-ascii');

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
