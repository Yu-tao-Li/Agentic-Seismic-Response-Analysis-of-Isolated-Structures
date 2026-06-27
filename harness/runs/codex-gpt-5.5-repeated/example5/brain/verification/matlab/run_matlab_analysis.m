function run_matlab_analysis()
% Nonlinear 4-DOF isolated shear-building analysis.
% Units: kg, m, s, N. Coordinates are relative to the ground.

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(scriptDir, '..', '..');
gmFile = fullfile(rootDir, 'input_data', 'Northridge_01_NO_968.txt');

dt = 0.01;
g0 = 9.80665;
targetPga = 0.40 * g0;

m = [2.50e5; 2.70e5; 2.70e5; 1.80e5];
k = [50e6; 245e6; 195e6; 98e6];
fy = [500e3; 1225e3; 975e3; 490e3];
b = [0.05; 0.0; 0.0; 0.0];
xi = 0.05;

raw = load(gmFile);
raw = raw(:);
rawPeak = max(abs(raw));
if rawPeak <= 0
    error('Ground motion file has zero peak acceleration.');
end
ag = raw * (targetPga / rawPeak);
n = numel(ag);
t = (0:n-1)' * dt;

M = diag(m);
K0 = assemble_story_matrix(k);
[modeShape, omega2] = eig(K0, M); %#ok<ASGLU>
omega = sqrt(sort(diag(omega2)));
omega = omega(isfinite(omega) & omega > 0);
if numel(omega) < 4
    error('Expected four positive elastic frequencies.');
end
A = [1/(2*omega(1)), omega(1)/2; 1/(2*omega(4)), omega(4)/2];
coeff = A \ [xi; xi];
alphaM = coeff(1);
betaK = coeff(2);
C = alphaM * M + betaK * K0;

gammaN = 0.5;
betaN = 0.25;
a0 = 1/(betaN*dt^2);
a2 = 1/(betaN*dt);
a3 = 1/(2*betaN) - 1;
cv = gammaN/(betaN*dt);

u = zeros(4, n);
v = zeros(4, n);
acc = zeros(4, n);
storyForce = zeros(4, n);
storyDef = zeros(4, n);
storyPlastic = zeros(4, n);
yielding = false(4, n);
iterCount = zeros(n, 1);
resNorm = zeros(n, 1);
converged = true(n, 1);

state = init_material_state(4);
[global0, ~, trialState, def0, local0] = evaluate_stories(u(:, 1), state, k, fy, b);
state = trialState;
storyForce(:, 1) = local0;
storyDef(:, 1) = def0;
storyPlastic(:, 1) = [state.plastic]';
acc(:, 1) = M \ (-M * ones(4, 1) * ag(1) - C * v(:, 1) - global0);

maxIter = 60;
failedSteps = 0;
for istep = 2:n
    p = -M * ones(4, 1) * ag(istep);
    uPrev = u(:, istep-1);
    vPrev = v(:, istep-1);
    aPrev = acc(:, istep-1);

    uTrial = uPrev + dt*vPrev + dt^2*(0.5 - betaN)*aPrev;
    stepConverged = false;
    lastRes = inf;
    trialState = state;
    fInt = zeros(4, 1);
    localForce = zeros(4, 1);
    ktStory = k;
    defTrial = zeros(4, 1);

    for iter = 1:maxIter
        [fInt, ktStory, trialState, defTrial, localForce] = evaluate_stories(uTrial, state, k, fy, b);
        aTrial = a0*(uTrial - uPrev) - a2*vPrev - a3*aPrev;
        vTrial = vPrev + dt*((1 - gammaN)*aPrev + gammaN*aTrial);
        residual = M*aTrial + C*vTrial + fInt - p;
        Ktan = assemble_story_matrix(ktStory);
        Keff = a0*M + cv*C + Ktan;
        du = -Keff \ residual;
        uTrial = uTrial + du;

        lastRes = norm(residual);
        forceTol = max(1.0, 1.0e-7 * max(1.0, norm(p)));
        dispTol = 1.0e-10 * max(1.0, norm(uTrial));
        if lastRes <= forceTol && norm(du) <= dispTol
            stepConverged = true;
            break;
        end
    end

    if ~stepConverged
        failedSteps = failedSteps + 1;
        converged(istep) = false;
    end

    [fInt, ktStory, trialState, defTrial, localForce] = evaluate_stories(uTrial, state, k, fy, b); %#ok<ASGLU>
    aTrial = a0*(uTrial - uPrev) - a2*vPrev - a3*aPrev;
    vTrial = vPrev + dt*((1 - gammaN)*aPrev + gammaN*aTrial);

    u(:, istep) = uTrial;
    v(:, istep) = vTrial;
    acc(:, istep) = aTrial;
    state = trialState;
    storyForce(:, istep) = localForce;
    storyDef(:, istep) = defTrial;
    storyPlastic(:, istep) = [state.plastic]';
    yielding(:, istep) = abs(storyPlastic(:, istep)) > 1.0e-12;
    iterCount(istep) = iter;
    resNorm(istep) = lastRes;
end

topDisp = u(4, :)';
topAbsAcc = acc(4, :)' + ag;
isoDisp = u(1, :)';
baseShear = storyForce(1, :)';
loopArea = trapz(isoDisp, baseShear);

writematrix([t, topDisp], fullfile(scriptDir, 'topStoDisIso2.txt'), 'Delimiter', 'tab');
writematrix([t, topAbsAcc], fullfile(scriptDir, 'topStoAccIso2.txt'), 'Delimiter', 'tab');
writematrix([t, isoDisp, baseShear], fullfile(scriptDir, 'isolation_hysteresis.txt'), 'Delimiter', 'tab');
writematrix([t, storyDef', storyForce'], fullfile(scriptDir, 'story_hysteresis_all.txt'), 'Delimiter', 'tab');

response = struct();
response.software = 'MATLAB';
response.generated_by = mfilename;
response.units = struct('length', 'm', 'mass', 'kg', 'time', 's', 'force', 'N', 'acceleration', 'm/s^2');
response.model = struct('masses_kg', m', 'stiffness_N_per_m', k', 'yield_N', fy', ...
    'post_yield_ratio', b', 'damping_ratio', xi, 'rayleigh_modes', [1, 4], ...
    'rayleigh_alphaM', alphaM, 'rayleigh_betaKinit', betaK);
response.input = struct('ground_motion_file', gmFile, 'dt', dt, 'sample_count', n, ...
    'raw_peak', rawPeak, 'target_pga_mps2', targetPga, ...
    'normalized_peak_mps2', max(abs(ag)));
response.time = t';
response.ground_acceleration = ag';
response.top_story_displacement = topDisp';
response.top_story_acceleration = topAbsAcc';
response.isolation_displacement = isoDisp';
response.base_shear = baseShear';
response.isolation_hysteresis = struct('displacement', isoDisp', 'force', baseShear', 'loop_area', loopArea);
response.story_displacement = storyDef';
response.story_force = storyForce';
response.story_plastic_deformation = storyPlastic';
response.yielding = yielding';
response.convergence = struct('failed_step_count', failedSteps, 'converged', converged', ...
    'iterations', iterCount', 'residual_norm', resNorm', 'max_iterations', max(iterCount));
response.peaks = struct('top_displacement_abs_max', max(abs(topDisp)), ...
    'top_acceleration_abs_max', max(abs(topAbsAcc)), ...
    'isolation_displacement_abs_max', max(abs(isoDisp)), ...
    'base_shear_abs_max', max(abs(baseShear)));

jsonText = jsonencode(response, PrettyPrint=true);
fid = fopen(fullfile(scriptDir, 'response.json'), 'w');
if fid < 0
    error('Could not open response.json for writing.');
end
fwrite(fid, jsonText, 'char');
fclose(fid);

fprintf('MATLAB nonlinear analysis complete. Samples: %d, failed steps: %d\n', n, failedSteps);
fprintf('Rayleigh alphaM = %.12g, betaKinit = %.12g\n', alphaM, betaK);
fprintf('Peaks: top displacement %.12g m, top acceleration %.12g m/s^2, isolation displacement %.12g m, base shear %.12g N\n', ...
    response.peaks.top_displacement_abs_max, response.peaks.top_acceleration_abs_max, ...
    response.peaks.isolation_displacement_abs_max, response.peaks.base_shear_abs_max);
end

function K = assemble_story_matrix(storyStiffness)
K = zeros(4, 4);
for i = 1:4
    e = zeros(4, 1);
    e(i) = 1;
    if i > 1
        e(i-1) = -1;
    end
    K = K + storyStiffness(i) * (e * e');
end
end

function state = init_material_state(n)
blank = struct('plastic', 0.0, 'backstress', 0.0);
state = repmat(blank, n, 1);
end

function [globalForce, ktStory, trialState, def, storyForce] = evaluate_stories(u, committedState, k, fy, b)
def = [u(1); u(2)-u(1); u(3)-u(2); u(4)-u(3)];
trialState = committedState;
storyForce = zeros(4, 1);
ktStory = zeros(4, 1);
for i = 1:4
    [storyForce(i), ktStory(i), trialState(i)] = material_trial(def(i), committedState(i), k(i), fy(i), b(i));
end
globalForce = [storyForce(1)-storyForce(2); storyForce(2)-storyForce(3); ...
    storyForce(3)-storyForce(4); storyForce(4)];
end

function [force, tangent, newState] = material_trial(deformation, oldState, k, fy, b)
if b >= 1.0
    error('Post-yield ratio must be less than 1.0.');
end
H = b * k / max(1.0e-16, 1.0 - b);
trialForce = k * (deformation - oldState.plastic);
xi = trialForce - oldState.backstress;
yieldValue = abs(xi) - fy;
newState = oldState;
if yieldValue <= 1.0e-10 * max(1.0, fy)
    force = trialForce;
    tangent = k;
    return;
end
direction = sign(xi);
if direction == 0
    direction = 1;
end
dgamma = yieldValue / (k + H);
newState.plastic = oldState.plastic + dgamma * direction;
newState.backstress = oldState.backstress + H * dgamma * direction;
force = trialForce - k * dgamma * direction;
tangent = k * H / (k + H);
end
