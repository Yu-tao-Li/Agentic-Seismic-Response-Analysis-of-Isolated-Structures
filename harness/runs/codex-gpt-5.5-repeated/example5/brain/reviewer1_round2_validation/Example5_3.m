function response = Example5_3(forceRtol, displacementRtol, stateRtol, writeOutputs)
% Corrected Example 5 solver with dimensionless numerical criteria.
% Units: kg, m, s, N. Structural responses are relative to the ground.

if nargin < 1 || isempty(forceRtol)
    forceRtol = 1.0e-7;
end
if nargin < 2 || isempty(displacementRtol)
    displacementRtol = 1.0e-10;
end
if nargin < 3 || isempty(stateRtol)
    stateRtol = 1.0e-10;
end
if nargin < 4
    writeOutputs = true;
end

scriptDir = fileparts(mfilename('fullpath'));
gmFile = fullfile(scriptDir, '..', 'input_data', 'Northridge_01_NO_968.txt');

dt = 0.01;
g0 = 9.81;
targetPga = 0.40 * g0;

raw = load(gmFile);
raw = raw(:, 1);
rawPeak = max(abs(raw));
if rawPeak <= 0
    error('Ground-motion record has zero peak acceleration.');
end
ag = raw * (targetPga / rawPeak);
time = (0:numel(ag)-1)' * dt;

m = [250; 270; 270; 180] * 1.0e3;
k = [50; 245; 195; 98] * 1.0e6;
fy = [500; 1225; 975; 490] * 1.0e3;
b = [0.05; 0.0; 0.0; 0.0];
nStories = numel(m);

M = diag(m);
K0 = computeStiffnessMatrix(k);
[~, D] = eig(K0, M);
omega = sort(sqrt(diag(D)));
dampingRatio = 0.05;
alphaM = 2 * dampingRatio * omega(1) * omega(2) / (omega(1) + omega(2));
betaK = 2 * dampingRatio / (omega(1) + omega(2));
C = alphaM * M + betaK * K0;

betaN = 0.25;
gammaN = 0.5;
influence = ones(nStories, 1);

u = zeros(nStories, 1);
v = zeros(nStories, 1);
restoringForce = zeros(nStories, 1);
a = M \ (-M * influence * ag(1) - C * v - restoringForce);

nSteps = numel(time);
topDisplacement = zeros(nSteps, 1);
topAcceleration = zeros(nSteps, 1);
isolationDisplacement = zeros(nSteps, 1);
baseShear = zeros(nSteps, 1);
convergenceFailures = false(nSteps, 1);
iterationCount = zeros(nSteps, 1);
residualRatioHistory = zeros(nSteps, 1);
correctionRatioHistory = zeros(nSteps, 1);
stateRatioHistory = zeros(nSteps, nStories);

topDisplacement(1) = u(end);
topAcceleration(1) = a(end);
isolationDisplacement(1) = u(1);
baseShear(1) = sum(restoringForce);

Keff = K0 + M / (betaN * dt^2) + gammaN * C / (betaN * dt);
maxIterations = 100;

for step = 1:nSteps-1
    deltaP = -M * influence * (ag(step+1) - ag(step));
    deltaPeff = deltaP ...
        + (M / (betaN * dt) + gammaN * C / betaN) * v ...
        + (M / (2 * betaN) + dt * (gammaN / (2 * betaN) - 1) * C) * a;

    uTrial = u;
    restoringTrial = restoringForce;
    residual = deltaPeff;
    converged = false;
    stateRatio = zeros(nStories, 1);

    for iteration = 1:maxIterations
        deltaU = Keff \ residual;
        uNext = uTrial + deltaU;

        [restoringNext, tangentStory, stateRatio] = computeRsStiffness( ...
            restoringTrial, uTrial, uNext, fy, k, b, stateRtol);
        Ktrial = computeStiffnessMatrix(tangentStory);
        KeffNext = Ktrial + M / (betaN * dt^2) + gammaN * C / (betaN * dt);

        deltaF = restoringNext - restoringTrial ...
            + gammaN * C * deltaU / (betaN * dt) ...
            + M * deltaU / (betaN * dt^2);
        residual = residual - deltaF;

        uTrial = uNext;
        restoringTrial = restoringNext;
        Keff = KeffNext;

        forceScale = max([norm(deltaPeff, 1), norm(restoringTrial, 1), norm(fy, 1), eps]);
        residualRatio = norm(residual, 1) / forceScale;
        displacementScale = max([norm(uTrial, inf), max(fy ./ k), eps]);
        correctionRatio = norm(deltaU, inf) / displacementScale;
        if residualRatio <= forceRtol && correctionRatio <= displacementRtol
            converged = true;
            break;
        end
    end

    iterationCount(step+1) = iteration;
    residualRatioHistory(step+1) = residualRatio;
    correctionRatioHistory(step+1) = correctionRatio;
    stateRatioHistory(step+1, :) = stateRatio';
    convergenceFailures(step+1) = ~converged;

    deltaU = uTrial - u;
    deltaV = gammaN * deltaU / (betaN * dt) - gammaN * v / betaN ...
        + (1 - gammaN / (2 * betaN)) * dt * a;
    deltaA = deltaU / (betaN * dt^2) - v / (betaN * dt) - a / (2 * betaN);

    u = uTrial;
    v = v + deltaV;
    a = a + deltaA;
    restoringForce = restoringTrial;

    topDisplacement(step+1) = u(end);
    topAcceleration(step+1) = a(end);
    isolationDisplacement(step+1) = u(1);
    baseShear(step+1) = sum(restoringForce);
end

response = struct();
response.software = 'MATLAB';
response.generated_by = mfilename;
response.units = struct('length', 'm', 'mass', 'kg', 'time', 's', ...
    'force', 'N', 'acceleration', 'm/s^2');
response.criteria = struct('force_residual_relative_tolerance', forceRtol, ...
    'displacement_correction_relative_tolerance', displacementRtol, ...
    'state_transition_relative_tolerance', stateRtol, ...
    'force_residual_definition', 'norm(residual,1)/max(norm(deltaPeff,1),norm(restoringForce,1),norm(fy,1))', ...
    'displacement_correction_definition', 'norm(deltaU,inf)/max(norm(uTrial,inf),max(fy./k))', ...
    'state_transition_definition', 'abs(shear-elasticTrial)/max(abs(shear),abs(elasticTrial),fy)');
response.model = struct('masses_kg', m', 'stiffness_N_per_m', k', ...
    'yield_N', fy', 'post_yield_ratio', b', 'damping_ratio', dampingRatio, ...
    'rayleigh_modes', [1, 2], 'rayleigh_alphaM', alphaM, 'rayleigh_betaKinit', betaK);
response.input = struct('ground_motion_file', gmFile, 'dt', dt, ...
    'sample_count', nSteps, 'target_pga_mps2', targetPga);
response.convergence = struct('failed_step_count', sum(convergenceFailures), ...
    'max_iterations', max(iterationCount), ...
    'max_final_residual_ratio', max(residualRatioHistory), ...
    'max_final_correction_ratio', max(correctionRatioHistory));
response.peaks = struct('top_displacement_abs_max', max(abs(topDisplacement)), ...
    'top_acceleration_abs_max', max(abs(topAcceleration)), ...
    'isolation_displacement_abs_max', max(abs(isolationDisplacement)), ...
    'base_shear_abs_max', max(abs(baseShear)), ...
    'hysteresis_loop_area_abs', abs(trapz(isolationDisplacement, baseShear)));
response.series = struct('time', time, 'top_displacement', topDisplacement, ...
    'top_acceleration', topAcceleration, ...
    'isolation_displacement', isolationDisplacement, 'base_shear', baseShear, ...
    'residual_ratio', residualRatioHistory, ...
    'correction_ratio', correctionRatioHistory, 'state_ratio', stateRatioHistory);

if writeOutputs
    writematrix([time, topDisplacement], fullfile(scriptDir, 'topStoDisIso2.txt'), 'Delimiter', 'tab');
    writematrix([time, topAcceleration], fullfile(scriptDir, 'topStoAccIso2.txt'), 'Delimiter', 'tab');
    writematrix([time, isolationDisplacement, baseShear], ...
        fullfile(scriptDir, 'hysteresis_layer1.txt'), 'Delimiter', 'tab');
    writematrix([time, double(convergenceFailures), iterationCount, ...
        residualRatioHistory, correctionRatioHistory], ...
        fullfile(scriptDir, 'convergence_log.txt'), 'Delimiter', 'tab');

    jsonResponse = rmfield(response, 'series');
    fid = fopen(fullfile(scriptDir, 'response_summary.json'), 'w');
    if fid < 0
        error('Could not write response_summary.json.');
    end
    fwrite(fid, jsonencode(jsonResponse, PrettyPrint=true), 'char');
    fclose(fid);
end

fprintf(['Example 5 complete: force rtol %.1e, displacement rtol %.1e, ' ...
    'state rtol %.1e, failed steps %d.\n'], forceRtol, displacementRtol, ...
    stateRtol, response.convergence.failed_step_count);
end
