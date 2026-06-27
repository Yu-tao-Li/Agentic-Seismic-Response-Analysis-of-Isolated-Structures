% Three-story linear shear building — Newmark-beta (average acceleration)
% Rayleigh damping xi=0.05, ground motion normalized to 0.40g
clear; close all;

nStories = 3;
m       = 1000;
k_story = 500e3;
xi      = 0.05;
dt      = 0.02;
g       = 9.81;
targetPGA = 0.40 * g;

%% Mass and stiffness matrices
M = m * eye(nStories);
K = zeros(nStories);
for s = 1:nStories
    if s == 1
        K(s,s) = k_story * 2;
        K(s,s+1) = -k_story;
    elseif s == nStories
        K(s,s-1) = -k_story;
        K(s,s) = k_story;
    else
        K(s,s-1) = -k_story;
        K(s,s)   = k_story * 2;
        K(s,s+1) = -k_story;
    end
end

%% Rayleigh damping from first two modes
[evec, ev] = eig(K, M);
omega = sort(sqrt(diag(ev)));
omega1 = omega(1);
omega2 = omega(2);
alpha = 2 * xi * omega1 * omega2 / (omega1 + omega2);
beta  = 2 * xi / (omega1 + omega2);
C = alpha * M + beta * K;

fprintf('omega1=%.6f  omega2=%.6f  omega3=%.6f\n', omega1, omega2, omega(3));
fprintf('alpha=%.6f  beta=%.8f\n', alpha, beta);

%% Load and normalize ground motion
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
gmFile = fullfile(scriptDir, '..', '..', 'input_data', 'GM1.txt');
rawGM = load(gmFile);
np = length(rawGM);
rawPGA = max(abs(rawGM));
scale = targetPGA / rawPGA;
ugddot = rawGM * scale;
fprintf('N=%d  rawPGA=%.6f  scale=%.6f  PGA=%.6f\n', np, rawPGA, scale, max(abs(ugddot)));

%% Newmark-beta average acceleration (gamma=0.5, beta=0.25)
gamma_nb = 0.5;
beta_nb  = 0.25;

a0 = 1 / (beta_nb * dt^2);
a1 = gamma_nb / (beta_nb * dt);
a2 = 1 / (beta_nb * dt);
a3 = 1/(2*beta_nb) - 1;
a4 = gamma_nb/beta_nb - 1;
a5 = dt/2 * (gamma_nb/beta_nb - 2);
a6 = dt * (1 - gamma_nb);
a7 = gamma_nb * dt;

Keff = K + a0*M + a1*C;
[L_eff, U_eff, P_eff] = lu(Keff);

infVec = ones(nStories, 1);

u   = zeros(nStories, 1);
ud  = zeros(nStories, 1);
udd = zeros(nStories, 1);

time   = zeros(np, 1);
topDis = zeros(np, 1);
topAcc = zeros(np, 1);

for i = 1:np
    time(i)   = (i-1) * dt;
    topDis(i) = u(nStories);
    topAcc(i) = udd(nStories);
    if i == np, break; end

    P_ext = -M * infVec * ugddot(i+1);
    F_eff = P_ext + M*(a0*u + a2*ud + a3*udd) + C*(a1*u + a4*ud + a5*udd);
    u_new = U_eff \ (L_eff \ (P_eff * F_eff));
    udd_new = a0*(u_new - u) - a2*ud - a3*udd;
    ud_new  = ud + a6*udd + a7*udd_new;

    u = u_new; ud = ud_new; udd = udd_new;
end

%% Write topStoDis.txt, topStoAcc.txt
writeMatrix(fullfile(pwd, 'topStoDis.txt'), [time, topDis], '%.10e\n');
writeMatrix(fullfile(pwd, 'topStoAcc.txt'), [time, topAcc], '%.10e\n');

fid = fopen(fullfile(pwd, 'response_output.json'), 'w');
fprintf(fid, '{\n');
fprintf(fid, '  "software": "MATLAB",\n');
fprintf(fid, '  "time_step": %.6f,\n', dt);
fprintf(fid, '  "n_samples": %d,\n', np);
fprintf(fid, '  "input_pga_normalized_m_s2": %.6f,\n', targetPGA);
fprintf(fid, '  "raw_pga_m_s2": %.6f,\n', rawPGA);
fprintf(fid, '  "scale_factor": %.6f,\n', scale);
fprintf(fid, '  "top_displacement_peak_m": %.6e,\n', max(abs(topDis)));
fprintf(fid, '  "top_acceleration_peak_m_s2": %.6e,\n', max(abs(topAcc)));
fprintf(fid, '  "omega1_rad_s": %.6f,\n', omega1);
fprintf(fid, '  "omega2_rad_s": %.6f,\n', omega2);
fprintf(fid, '  "omega3_rad_s": %.6f,\n', omega(3));
fprintf(fid, '  "rayleigh_alpha": %.6f,\n', alpha);
fprintf(fid, '  "rayleigh_beta": %.8f\n', beta);
fprintf(fid, '}\n');
fclose(fid);

fprintf('Peak top displacement: %.6e m\n', max(abs(topDis)));
fprintf('Peak top acceleration: %.6e m/s^2\n', max(abs(topAcc)));

function writeMatrix(fname, data, fmt)
    fid = fopen(fname, 'w');
    for row = 1:size(data,1)
        for col = 1:size(data,2)
            fprintf(fid, fmt, data(row,col));
        end
    end
    fclose(fid);
end
