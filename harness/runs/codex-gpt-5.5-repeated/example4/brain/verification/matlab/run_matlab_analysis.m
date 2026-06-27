function run_matlab_analysis()
% Nonlinear four-DOF isolated shear building verification model.
% Outputs are generated in this folder and are intended for cross-validation.

clearvars;
clc;

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fileparts(scriptDir));
inputFile = fullfile(rootDir, 'input_data', 'Northridge_01_NO_968.txt');
logFile = fullfile(scriptDir, 'matlab_internal_log.txt');
diary(logFile);
diary on;

try
    fprintf('MATLAB nonlinear MDOF analysis started: %s\n', datestr(now, 31));
    raw = load(inputFile);
    raw = raw(:);
    if isempty(raw)
        error('Input ground motion file is empty: %s', inputFile);
    end

    dt = 0.01;
    g = 9.81;
    targetPGA = 0.40 * g;
    scale = targetPGA / max(abs(raw));
    ag = raw * scale;
    n = numel(ag);
    time = (0:n-1)' * dt;

    model = model_parameters();
    [alphaM, betaK, omegas] = rayleigh_from_initial_modes(model.M, model.K0, 0.05);
    C = alphaM * model.M + betaK * model.K0;
    C(1, 1) = C(1, 1) + model.cIso;

    newmark.gamma = 0.5;
    newmark.beta = 0.25;
    tol = 1.0e-7;
    maxIter = 50;

    u = zeros(4, n);
    v = zeros(4, n);
    a = zeros(4, n);
    storyForce = zeros(4, n);
    baseShear = zeros(n, 1);
    tangentCode = zeros(4, n);
    yieldedInstant = false(n, 3);
    yieldedCumulative = false(n, 3);
    firstYieldTime = nan(3, 1);
    iterCount = zeros(n, 1);
    converged = true(n, 1);
    failedSteps = 0;

    plastic = zeros(3, 1);
    [fint0, q0, kt0, instYield0, plastic0] = restoring_force(u(:, 1), plastic, model);
    plastic = plastic0;
    storyForce(:, 1) = q0;
    tangentCode(:, 1) = kt0 > 0;
    yieldedInstant(1, :) = instYield0(:)';
    yieldedCumulative(1, :) = instYield0(:)';
    for s = 1:3
        if instYield0(s)
            firstYieldTime(s) = time(1);
        end
    end
    a(:, 1) = model.M \ (-model.M * ones(4, 1) * ag(1) - C * v(:, 1) - fint0);
    baseShear(1) = q0(1);

    for i = 2:n
        pNext = -model.M * ones(4, 1) * ag(i);
        [u(:, i), v(:, i), a(:, i), plastic, q, kt, instYield, ok, iters] = ...
            newmark_step(u(:, i-1), v(:, i-1), a(:, i-1), plastic, pNext, model, C, newmark, dt, tol, maxIter);

        converged(i) = ok;
        iterCount(i) = iters;
        if ~ok
            failedSteps = failedSteps + 1;
        end
        storyForce(:, i) = q;
        baseShear(i) = q(1);
        tangentCode(:, i) = kt > 0;
        yieldedInstant(i, :) = instYield(:)';
        yieldedCumulative(i, :) = yieldedCumulative(i-1, :) | instYield(:)';
        for s = 1:3
            if instYield(s) && isnan(firstYieldTime(s))
                firstYieldTime(s) = time(i);
            end
        end
    end

    topDisp = u(4, :)';
    isoDisp = u(1, :)';
    topRelAcc = a(4, :)';
    topAbsAcc = topRelAcc + ag;

    writematrix([time, topDisp], fullfile(scriptDir, 'topStoDisIso1.txt'), 'Delimiter', 'tab');
    writematrix([time, topAbsAcc], fullfile(scriptDir, 'topStoAccIso1.txt'), 'Delimiter', 'tab');
    writematrix([time, isoDisp], fullfile(scriptDir, 'isoDisIso1.txt'), 'Delimiter', 'tab');
    writematrix([time, baseShear], fullfile(scriptDir, 'baseShearIso1.txt'), 'Delimiter', 'tab');

    response = struct();
    response.software = 'MATLAB';
    response.generated_at = char(datetime('now', 'TimeZone', 'local', 'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
    response.input = struct('file', inputFile, 'raw_count', n, 'dt', dt, ...
        'scale_to_mps2', scale, 'normalized_pga_mps2', max(abs(ag)), 'target_pga_g', 0.40);
    response.model = export_model(model, alphaM, betaK, omegas);
    response.method = struct('integrator', 'Newmark average acceleration', ...
        'nonlinear_update', 'Newton iteration with elastic-perfectly-plastic return mapping', ...
        'tolerance', tol, 'max_iterations', maxIter);
    response.time = time';
    response.top_story_displacement = topDisp';
    response.top_story_acceleration = topAbsAcc';
    response.top_story_relative_acceleration = topRelAcc';
    response.isolation_displacement = isoDisp';
    response.restoring_force = storyForce';
    response.base_shear = baseShear';
    response.yielding = struct('instant', double(yieldedInstant), ...
        'cumulative', double(yieldedCumulative), ...
        'first_yield_time', firstYieldTime');
    response.convergence = struct('failed_step_count', failedSteps, ...
        'converged', double(converged'), 'iterations', iterCount');
    response.summary = make_summary(time, topDisp, topAbsAcc, isoDisp, baseShear, firstYieldTime, failedSteps);

    write_json(fullfile(scriptDir, 'response_matlab.json'), response);
    fprintf('MATLAB nonlinear MDOF analysis completed. Failed steps: %d\n', failedSteps);
catch ME
    fprintf(2, 'MATLAB analysis failed: %s\n', ME.message);
    for k = 1:numel(ME.stack)
        fprintf(2, '  at %s line %d\n', ME.stack(k).name, ME.stack(k).line);
    end
    failure = struct('software', 'MATLAB', 'status', 'failed', 'message', ME.message);
    write_json(fullfile(scriptDir, 'response_matlab.json'), failure);
    diary off;
    rethrow(ME);
end

diary off;
end

function model = model_parameters()
model.m = [2.50e5; 2.70e5; 2.70e5; 1.80e5];
model.k = [50e6; 245e6; 195e6; 98e6];
model.fy = [Inf; 1225e3; 975e3; 490e3];
model.b = [NaN; 0.0; 0.0; 0.0];
model.cIso = 1000e3;
model.M = diag(model.m);
B = [1 0 0 0; -1 1 0 0; 0 -1 1 0; 0 0 -1 1];
model.B = B;
model.K0 = B' * diag(model.k) * B;
end

function [alphaM, betaK, omegas] = rayleigh_from_initial_modes(M, K, zeta)
[~, D] = eig(K, M);
omegas = sqrt(sort(diag(D)));
w1 = omegas(1);
w2 = omegas(end);
A = [1/(2*w1), w1/2; 1/(2*w2), w2/2];
x = A \ [zeta; zeta];
alphaM = x(1);
betaK = x(2);
end

function [uNext, vNext, aNext, plasticNext, qNext, ktNext, instYield, ok, iters] = ...
    newmark_step(uPrev, vPrev, aPrev, plasticPrev, pNext, model, C, nm, dt, tol, maxIter)

a0 = 1 / (nm.beta * dt^2);
a1 = nm.gamma / (nm.beta * dt);
a2 = 1 / (nm.beta * dt);
a3 = 1 / (2 * nm.beta) - 1;
a4 = nm.gamma / nm.beta - 1;
a5 = dt * (nm.gamma / (2 * nm.beta) - 1);

uTrial = uPrev + dt * vPrev + dt^2 * (0.5 - nm.beta) * aPrev;
ok = false;
plasticTrial = plasticPrev;
qNext = zeros(4, 1);
ktNext = zeros(4, 1);
instYield = false(3, 1);

for iters = 1:maxIter
    aTrial = a0 * (uTrial - uPrev) - a2 * vPrev - a3 * aPrev;
    vTrial = vPrev + dt * ((1 - nm.gamma) * aPrev + nm.gamma * aTrial);
    [fint, q, kt, instYieldTrial, plasticCandidate] = restoring_force(uTrial, plasticPrev, model);
    residual = model.M * aTrial + C * vTrial + fint - pNext;
    keff = model.B' * diag(kt) * model.B + a0 * model.M + a1 * C;
    du = -keff \ residual;
    uTrial = uTrial + du;
    forceScale = max([norm(pNext, inf), norm(fint, inf), 1.0]);
    if norm(residual, inf) / forceScale < tol && norm(du, inf) < 1.0e-9
        ok = true;
        plasticTrial = plasticCandidate;
        qNext = q;
        ktNext = kt;
        instYield = instYieldTrial;
        break;
    end
    plasticTrial = plasticCandidate;
    qNext = q;
    ktNext = kt;
    instYield = instYieldTrial;
end

uNext = uTrial;
aNext = a0 * (uNext - uPrev) - a2 * vPrev - a3 * aPrev;
vNext = vPrev + dt * ((1 - nm.gamma) * aPrev + nm.gamma * aNext);
if ~ok
    [~, qNext, ktNext, instYield, plasticTrial] = restoring_force(uNext, plasticPrev, model);
end
plasticNext = plasticTrial;
end

function [fint, q, kt, instYield, plasticTrial] = restoring_force(u, plasticCommitted, model)
d = model.B * u;
q = zeros(4, 1);
kt = zeros(4, 1);
plasticTrial = plasticCommitted;
instYield = false(3, 1);

q(1) = model.k(1) * d(1);
kt(1) = model.k(1);

for j = 2:4
    idx = j - 1;
    fTrial = model.k(j) * (d(j) - plasticCommitted(idx));
    if fTrial > model.fy(j)
        q(j) = model.fy(j);
        kt(j) = model.b(j) * model.k(j);
        plasticTrial(idx) = d(j) - q(j) / model.k(j);
        instYield(idx) = true;
    elseif fTrial < -model.fy(j)
        q(j) = -model.fy(j);
        kt(j) = model.b(j) * model.k(j);
        plasticTrial(idx) = d(j) - q(j) / model.k(j);
        instYield(idx) = true;
    else
        q(j) = fTrial;
        kt(j) = model.k(j);
        instYield(idx) = abs(q(j)) >= 0.999999 * model.fy(j);
    end
end

fint = model.B' * q;
end

function summary = make_summary(time, topDisp, topAbsAcc, isoDisp, baseShear, firstYieldTime, failedSteps)
[peakTopDisp, i1] = max(abs(topDisp));
[peakTopAcc, i2] = max(abs(topAbsAcc));
[peakIsoDisp, i3] = max(abs(isoDisp));
[peakBaseShear, i4] = max(abs(baseShear));
summary = struct();
summary.peak_top_displacement_m = peakTopDisp;
summary.peak_top_displacement_time_s = time(i1);
summary.peak_top_acceleration_mps2 = peakTopAcc;
summary.peak_top_acceleration_g = peakTopAcc / 9.81;
summary.peak_top_acceleration_time_s = time(i2);
summary.peak_isolation_displacement_m = peakIsoDisp;
summary.peak_isolation_displacement_time_s = time(i3);
summary.peak_base_shear_N = peakBaseShear;
summary.peak_base_shear_time_s = time(i4);
summary.first_yield_time_s = firstYieldTime(:)';
summary.failed_step_count = failedSteps;
end

function out = export_model(model, alphaM, betaK, omegas)
out = struct();
out.masses_kg = model.m';
out.story_stiffness_N_per_m = model.k';
out.yield_forces_N = model.fy';
out.post_yield_ratios = model.b';
out.isolation_dashpot_Ns_per_m = model.cIso;
out.rayleigh = struct('zeta', 0.05, 'mode_pair', [1, 4], ...
    'alpha_mass', alphaM, 'beta_initial_stiffness', betaK, ...
    'c_matrix_Ns_per_m', alphaM * model.M + betaK * model.K0 + diag([model.cIso, 0, 0, 0]), ...
    'initial_omega_rad_per_s', omegas');
end

function write_json(path, data)
txt = jsonencode(data, 'PrettyPrint', true);
fid = fopen(path, 'w');
if fid < 0
    error('Cannot open JSON output for writing: %s', path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', txt);
end
