% MATLAB reference benchmark for Example 6.
% Three-story, two-bay elastic steel frame with Bouc-Wen base isolators.
% Writes response_output.json in this directory.

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end
outputPath = fullfile(scriptDir, 'response_output.json');

try
    payload = run_example6_matlab(scriptDir);
catch ME
    payload = unavailable_payload('failed', sprintf('%s: %s', ME.identifier, ME.message));
end

write_json(outputPath, payload);

function payload = run_example6_matlab(scriptDir)
    runRoot = fileparts(fileparts(scriptDir));
    inputPath = fullfile(runRoot, 'input_data', 'Northridge_01_NO_968.txt');
    historyCsv = fullfile(scriptDir, 'time_histories.csv');
    hysteresisCsv = fullfile(scriptDir, 'isolator_hysteresis.csv');
    normalizedPath = fullfile(scriptDir, 'normalized_ground_motion.txt');

    p = example6_params();
    if ~exist(inputPath, 'file')
        payload = unavailable_payload('missing_input', ['Required ground-motion file is missing: ', inputPath]);
        return;
    end

    raw = load(inputPath);
    raw = raw(:);
    if isempty(raw)
        payload = unavailable_payload('empty_input', ['Ground-motion file is empty: ', inputPath]);
        return;
    end
    rawPeak = max(abs(raw));
    if rawPeak <= 0.0
        payload = unavailable_payload('invalid_input', 'Ground-motion peak acceleration is zero.');
        return;
    end

    targetPga = p.target_pga_g * p.g;
    scale = targetPga / rawPeak;
    ag = raw * scale;
    nSteps = numel(ag);
    time = (0:nSteps-1)' * p.dt;
    write_vector(normalizedPath, ag);

    [Klinear, Kinitial, M, fullToReduced, isoIdx, floorIdx] = build_reduced_matrices(p);
    [fundamentalPeriod, modalPeriods, alphaM, betaK] = modal_and_damping(Kinitial, M, p.damping_ratio);
    C = alphaM * M + betaK * Kinitial;

    nd = size(M, 1);
    influence = zeros(nd, 1);
    massDiag = diag(M);
    influence(massDiag > 0.0) = 1.0;

    beta = p.newmark_beta;
    gamma = p.newmark_gamma;
    dt = p.dt;
    a0 = 1.0 / (beta * dt * dt);
    a1 = gamma / (beta * dt);
    tol = 1.0e-8;
    maxIter = 50;

    q = zeros(nd, 1);
    v = zeros(nd, 1);
    acc = zeros(nd, 1);
    isoState.u = zeros(3, 1);
    isoState.z = zeros(3, 1);

    roof = zeros(nSteps, 1);
    story = zeros(nSteps, 3);
    isoMean = zeros(nSteps, 1);
    baseShear = zeros(nSteps, 1);
    isoDispHist = zeros(nSteps, 3);
    isoForceHist = zeros(nSteps, 3);
    isoZHist = zeros(nSteps, 3);
    iterations = zeros(nSteps, 1);
    residualRatios = zeros(nSteps, 1);
    failedSteps = [];

    [resp, isoForce0] = collect_response(q, isoState, isoIdx, floorIdx, p);
    roof(1) = resp.roof;
    story(1, :) = resp.storyDrifts;
    isoMean(1) = resp.isolation;
    baseShear(1) = resp.baseShear;
    isoDispHist(1, :) = q(isoIdx)';
    isoForceHist(1, :) = isoForce0';
    isoZHist(1, :) = isoState.z';

    for istep = 2:nSteps
        qPred = q + dt * v + dt * dt * (0.5 - beta) * acc;
        vPred = v + dt * (1.0 - gamma) * acc;
        pExt = -M * influence * ag(istep);
        qTrial = q;
        converged = false;
        lastRatio = Inf;
        trialState = isoState;
        isoForce = zeros(3, 1);

        for iter = 1:maxIter
            accTrial = a0 * (qTrial - qPred);
            vTrial = vPred + gamma * dt * accTrial;
            [isoForce, isoTangent, trialState] = isolator_trial(qTrial(isoIdx), isoState, p);

            fint = Klinear * qTrial;
            kt = Klinear;
            for j = 1:3
                fint(isoIdx(j)) = fint(isoIdx(j)) + isoForce(j);
                kt(isoIdx(j), isoIdx(j)) = kt(isoIdx(j), isoIdx(j)) + isoTangent(j);
            end

            residual = M * accTrial + C * vTrial + fint - pExt;
            denom = max([1.0, norm(pExt, Inf), norm(fint, Inf)]);
            lastRatio = norm(residual, Inf) / denom;
            if lastRatio < tol
                converged = true;
                break;
            end

            keff = kt + a0 * M + a1 * C;
            dq = -keff \ residual;
            qTrial = qTrial + dq;
        end

        iterations(istep) = iter;
        residualRatios(istep) = lastRatio;
        if ~converged
            failedSteps(end + 1) = istep - 1; %#ok<AGROW>
            roof = roof(1:istep-1);
            story = story(1:istep-1, :);
            isoMean = isoMean(1:istep-1);
            baseShear = baseShear(1:istep-1);
            isoDispHist = isoDispHist(1:istep-1, :);
            isoForceHist = isoForceHist(1:istep-1, :);
            isoZHist = isoZHist(1:istep-1, :);
            iterations = iterations(1:istep-1);
            residualRatios = residualRatios(1:istep-1);
            time = time(1:istep-1);
            ag = ag(1:istep-1);
            break;
        end

        q = qTrial;
        acc = a0 * (q - qPred);
        v = vPred + gamma * dt * acc;
        isoState = trialState;

        [resp, isoForce] = collect_response(q, isoState, isoIdx, floorIdx, p);
        roof(istep) = resp.roof;
        story(istep, :) = resp.storyDrifts;
        isoMean(istep) = resp.isolation;
        baseShear(istep) = resp.baseShear;
        isoDispHist(istep, :) = q(isoIdx)';
        isoForceHist(istep, :) = isoForce';
        isoZHist(istep, :) = isoState.z';
    end

    write_history_csv(historyCsv, time, ag, roof, story, isoMean, baseShear);
    write_hysteresis_csv(hysteresisCsv, time, isoDispHist, isoForceHist, isoZHist);
    peaks = peak_quantities(roof, story, isoMean, baseShear, isoDispHist, isoForceHist);

    failedCount = numel(failedSteps);
    attempted = max(0, nSteps - 1);

    inputMeta.dt_s = p.dt;
    inputMeta.sample_count = numel(raw);
    inputMeta.raw_peak_abs = rawPeak;
    inputMeta.target_pga_m_per_s2 = targetPga;
    inputMeta.target_pga_g = p.target_pga_g;
    inputMeta.normalization_scale = scale;
    inputMeta.normalized_peak_abs_m_per_s2 = max(abs(raw * scale));
    inputMeta.source_file = 'input_data/Northridge_01_NO_968.txt';

    damping.type = 'Rayleigh';
    damping.target_modes = [1, 2];
    damping.target_damping_ratio = p.damping_ratio;
    damping.alpha_mass = alphaM;
    damping.beta_stiffness_initial = betaK;
    damping.period_source = 'MATLAB static condensation of initial elastic tangent stiffness';

    convergence.analysis = 'Newmark beta=0.25 gamma=0.5 with Newton-Raphson isolator iterations';
    convergence.test = 'relative infinity-norm residual';
    convergence.relative_tolerance = tol;
    convergence.max_iterations_per_step = maxIter;
    convergence.attempted_steps = attempted;
    convergence.completed_steps = numel(time) - 1;
    convergence.failed_step_indices = failedSteps;
    convergence.failed_step_count = failedCount;
    convergence.failed_step_rate = failedCount / max(1, attempted);
    convergence.iterations_per_completed_step = iterations(:)';
    convergence.relative_residual_per_completed_step = residualRatios(:)';

    payload.schema_version = 'example6-response-v1';
    payload.software = 'MATLAB';
    if failedCount == 0
        payload.status = 'completed';
    else
        payload.status = 'failed';
    end
    payload.generated_at_utc = utc_now_string();
    payload.units.length = 'm';
    payload.units.force = 'N';
    payload.units.mass = 'kg';
    payload.units.time = 's';
    payload.input = inputMeta;
    payload.model = p;
    payload.damping = damping;
    payload.fundamental_period_s = fundamentalPeriod;
    payload.modal_periods_s = modalPeriods(:)';
    payload.time_s = time(:)';
    payload.ground_accel_m_per_s2 = ag(:)';
    payload.roof_displacement_m = roof(:)';
    payload.story_drift_ratios.story_1 = story(:, 1)';
    payload.story_drift_ratios.story_2 = story(:, 2)';
    payload.story_drift_ratios.story_3 = story(:, 3)';
    payload.isolation_displacement_m = isoMean(:)';
    payload.total_base_shear_N = baseShear(:)';
    for j = 1:3
        name = sprintf('isolator_%d', j);
        payload.isolator_hysteresis.(name).deformation_m = isoDispHist(:, j)';
        payload.isolator_hysteresis.(name).force_N = isoForceHist(:, j)';
        payload.isolator_hysteresis.(name).boucwen_z_m = isoZHist(:, j)';
        payload.isolator_hysteresis.(name).loop_area_signed_Nm = peaks.isolator_loop_area_signed_Nm(j);
        payload.isolator_hysteresis.(name).loop_area_abs_Nm = peaks.isolator_loop_area_abs_Nm(j);
    end
    payload.peaks = peaks;
    payload.convergence = convergence;
    payload.assumptions = { ...
        'The 2D beam-column frame is assembled with ux, uy, and rz nodal DOFs before applying floor horizontal diaphragm transformations.', ...
        'Floor horizontal diaphragm constraints tie only ux at each floor; vertical translations and rotations remain frame DOFs.', ...
        'Floor mass is assigned only to horizontal translational DOFs, one third of each floor mass at each floor node.', ...
        'Horizontal isolators use F = alpha*k0*u + (1-alpha)*k0*z with a backward-Euler Bouc-Wen update for z.', ...
        'The Bouc-Wen update follows the OpenSees form z_{n+1} - z_n - du*(Ao - abs(z_{n+1})^n*(gamma + beta*sign(du*z_{n+1}))) = 0, with du = u_{n+1}-u_n, by local Newton iteration.', ...
        'Story-1 drift uses the average top-of-isolator displacement as the lower reference; isolation displacement is the mean of the three isolator deformations.', ...
        'Gravity load and P-Delta effects are intentionally excluded.' ...
    };
    payload.output_files.json = 'verification/matlab/response_output.json';
    payload.output_files.time_histories_csv = 'verification/matlab/time_histories.csv';
    payload.output_files.isolator_hysteresis_csv = 'verification/matlab/isolator_hysteresis.csv';
    payload.output_files.normalized_ground_motion = 'verification/matlab/normalized_ground_motion.txt';
end

function p = example6_params()
    p.g = 9.81;
    p.target_pga_g = 0.40;
    p.dt = 0.01;
    p.stories = 3;
    p.bays = 2;
    p.bay_width_m = 6.0;
    p.story_height_m = 3.6;
    p.E_Pa = 2.06e11;
    p.nu = 0.30;
    p.rho_kg_per_m3 = 7850.0;
    p.column_area_m2 = 0.020;
    p.column_I_m4 = 8.0e-4;
    p.beam_area_m2 = 0.015;
    p.beam_I_m4 = 4.5e-4;
    p.floor_mass_kg = 2.5e5;
    p.isolator_model = 'BoucWen';
    p.isolator_k0_N_per_m = 2.0e6;
    p.isolator_characteristic_yield_force_N = 1.0e4;
    p.isolator_characteristic_yield_displacement_m = 0.005;
    p.boucwen_alpha = 0.05;
    p.boucwen_n = 2.0;
    p.boucwen_beta = 2.0e4;
    p.boucwen_gamma = 2.0e4;
    p.boucwen_Ao = 1.0;
    p.boucwen_deltaA = 0.0;
    p.boucwen_deltaNu = 0.0;
    p.boucwen_deltaEta = 0.0;
    p.isolator_kv_N_per_m = 1.0e10;
    p.damping_ratio = 0.05;
    p.newmark_beta = 0.25;
    p.newmark_gamma = 0.50;
end

function [Klinear, Kinitial, M, fullToReduced, isoIdx, floorIdx] = build_reduced_matrices(p)
    nCols = 3;
    nLevels = 4;
    nNodes = nCols * nLevels;
    nFull = nNodes * 3;
    KlinearFull = zeros(nFull, nFull);
    KinitialFull = zeros(nFull, nFull);
    Mfull = zeros(nFull, nFull);

    coords = zeros(nNodes, 2);
    for level = 0:nLevels-1
        for col = 0:nCols-1
            nid = node_id(level, col);
            coords(nid, :) = [col * p.bay_width_m, level * p.story_height_m];
        end
    end

    for col = 0:nCols-1
        for level = 0:2
            ni = node_id(level, col);
            nj = node_id(level + 1, col);
            ke = frame_element_stiffness(p.E_Pa, p.column_area_m2, p.column_I_m4, coords(ni, :), coords(nj, :));
            ed = [node_dofs(ni), node_dofs(nj)];
            KlinearFull(ed, ed) = KlinearFull(ed, ed) + ke;
            KinitialFull(ed, ed) = KinitialFull(ed, ed) + ke;
        end
    end

    for level = 1:3
        for col = 0:1
            ni = node_id(level, col);
            nj = node_id(level, col + 1);
            ke = frame_element_stiffness(p.E_Pa, p.beam_area_m2, p.beam_I_m4, coords(ni, :), coords(nj, :));
            ed = [node_dofs(ni), node_dofs(nj)];
            KlinearFull(ed, ed) = KlinearFull(ed, ed) + ke;
            KinitialFull(ed, ed) = KinitialFull(ed, ed) + ke;
        end
    end

    for col = 0:nCols-1
        ux = dof_id(0, col, 1);
        uy = dof_id(0, col, 2);
        KinitialFull(ux, ux) = KinitialFull(ux, ux) + p.isolator_k0_N_per_m;
        KlinearFull(uy, uy) = KlinearFull(uy, uy) + p.isolator_kv_N_per_m;
        KinitialFull(uy, uy) = KinitialFull(uy, uy) + p.isolator_kv_N_per_m;
    end

    for level = 1:3
        for col = 0:nCols-1
            ux = dof_id(level, col, 1);
            Mfull(ux, ux) = Mfull(ux, ux) + p.floor_mass_kg / 3.0;
        end
    end

    keyMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    fullToReduced = zeros(nFull, 1);
    nextId = 0;
    for level = 0:nLevels-1
        for col = 0:nCols-1
            for comp = 1:3
                gd = dof_id(level, col, comp);
                if comp == 1 && level >= 1
                    key = sprintf('floor_ux_%d', level);
                else
                    key = sprintf('dof_%d', gd);
                end
                if ~isKey(keyMap, key)
                    nextId = nextId + 1;
                    keyMap(key) = nextId;
                end
                fullToReduced(gd) = keyMap(key);
            end
        end
    end

    T = zeros(nFull, nextId);
    for gd = 1:nFull
        T(gd, fullToReduced(gd)) = 1.0;
    end
    Klinear = T' * KlinearFull * T;
    Kinitial = T' * KinitialFull * T;
    M = T' * Mfull * T;

    isoIdx = zeros(3, 1);
    for col = 0:2
        isoIdx(col + 1) = fullToReduced(dof_id(0, col, 1));
    end
    floorIdx = zeros(3, 1);
    for level = 1:3
        floorIdx(level) = fullToReduced(dof_id(level, 0, 1));
    end
end

function [fundamentalPeriod, modalPeriods, alphaM, betaK] = modal_and_damping(Kinitial, M, zeta)
    dyn = find(diag(M) > 1.0e-9);
    all = (1:size(M, 1))';
    stat = setdiff(all, dyn);
    Kdd = Kinitial(dyn, dyn);
    Mdd = M(dyn, dyn);
    if isempty(stat)
        Kcond = Kdd;
    else
        Kds = Kinitial(dyn, stat);
        Ksd = Kinitial(stat, dyn);
        Kss = Kinitial(stat, stat);
        Kcond = Kdd - Kds * (Kss \ Ksd);
    end
    lambda = eig(Kcond, Mdd);
    lambda = sort(real(lambda(lambda > 1.0e-8)));
    omega = sqrt(lambda);
    modalPeriods = 2.0 * pi ./ omega;
    fundamentalPeriod = modalPeriods(1);
    alphaM = 2.0 * zeta * omega(1) * omega(2) / (omega(1) + omega(2));
    betaK = 2.0 * zeta / (omega(1) + omega(2));
end

function [force, tangent, trialState] = isolator_trial(disp, committedState, p)
    k0 = p.isolator_k0_N_per_m;
    alpha = p.boucwen_alpha;
    force = zeros(3, 1);
    tangent = zeros(3, 1);
    trialState.z = committedState.z;
    for j = 1:3
        [zNew, dzdu] = boucwen_update(disp(j), committedState.u(j), committedState.z(j), p);
        trialState.u(j) = disp(j);
        trialState.z(j) = zNew;
        force(j) = alpha * k0 * disp(j) + (1.0 - alpha) * k0 * zNew;
        tangent(j) = alpha * k0 + (1.0 - alpha) * k0 * dzdu;
    end
end

function [zNew, dzdu] = boucwen_update(uNew, uOld, zOld, p)
    n = p.boucwen_n;
    betaBW = p.boucwen_beta;
    gammaBW = p.boucwen_gamma;
    Ao = p.boucwen_Ao;
    du = uNew - uOld;
    z = zOld;
    tolLocal = 1.0e-12;
    for iterLocal = 1:30
        absz = abs(z);
        if absz < 1.0e-14
            absz = 1.0e-14;
        end
        psi = gammaBW + betaBW * sign_nonzero(du * z);
        phi = Ao - absz^n * psi;
        g = z - zOld - du * phi;
        dphiDz = -n * absz^(n - 1.0) * sign_nonzero(z) * psi;
        dgDz = 1.0 - du * dphiDz;
        dz = -g / dgDz;
        z = z + dz;
        if abs(dz) <= tolLocal * max(1.0, abs(z))
            break;
        end
    end
    absz = abs(z);
    if absz < 1.0e-14
        absz = 1.0e-14;
    end
    psi = gammaBW + betaBW * sign_nonzero(du * z);
    phi = Ao - absz^n * psi;
    dphiDz = -n * absz^(n - 1.0) * sign_nonzero(z) * psi;
    dgDz = 1.0 - du * dphiDz;
    dgDu = -phi;
    dzdu = -dgDu / dgDz;
    zNew = z;
end

function s = sign_nonzero(value)
    if value > 0.0
        s = 1.0;
    else
        s = -1.0;
    end
end

function [resp, isoForce] = collect_response(q, isoState, isoIdx, floorIdx, p)
    [isoForce, ~, ~] = isolator_trial(q(isoIdx), isoState, p);
    floorDisp = q(floorIdx);
    isoDisp = q(isoIdx);
    baseMean = mean(isoDisp);
    h = p.story_height_m;
    resp.roof = floorDisp(3);
    resp.storyDrifts = [(floorDisp(1) - baseMean) / h, (floorDisp(2) - floorDisp(1)) / h, (floorDisp(3) - floorDisp(2)) / h];
    resp.isolation = baseMean;
    resp.baseShear = sum(isoForce);
end

function peaks = peak_quantities(roof, story, isoMean, baseShear, isoDisp, isoForce)
    signed = zeros(1, 3);
    absArea = zeros(1, 3);
    for j = 1:3
        [signed(j), absArea(j)] = loop_area(isoForce(:, j), isoDisp(:, j));
    end
    peaks.peak_abs_roof_displacement_m = max(abs(roof));
    peaks.max_abs_story_drift_ratio = max(abs(story(:)));
    peaks.peak_abs_isolation_displacement_m = max(abs(isoMean));
    peaks.peak_abs_total_base_shear_N = max(abs(baseShear));
    peaks.isolator_loop_area_signed_Nm = signed;
    peaks.isolator_loop_area_abs_Nm = absArea;
    peaks.total_isolator_loop_area_signed_Nm = sum(signed);
    peaks.total_isolator_loop_area_abs_Nm = sum(absArea);
end

function [signedArea, absArea] = loop_area(force, disp)
    if numel(force) < 2
        signedArea = 0.0;
        absArea = 0.0;
        return;
    end
    inc = 0.5 * (force(2:end) + force(1:end-1)) .* diff(disp);
    signedArea = sum(inc);
    absArea = sum(abs(inc));
end

function ke = frame_element_stiffness(E, A, I, ci, cj)
    xi = ci(1); yi = ci(2);
    xj = cj(1); yj = cj(2);
    L = hypot(xj - xi, yj - yi);
    c = (xj - xi) / L;
    s = (yj - yi) / L;
    kl = [ ...
        A*E/L, 0, 0, -A*E/L, 0, 0; ...
        0, 12*E*I/L^3, 6*E*I/L^2, 0, -12*E*I/L^3, 6*E*I/L^2; ...
        0, 6*E*I/L^2, 4*E*I/L, 0, -6*E*I/L^2, 2*E*I/L; ...
        -A*E/L, 0, 0, A*E/L, 0, 0; ...
        0, -12*E*I/L^3, -6*E*I/L^2, 0, 12*E*I/L^3, -6*E*I/L^2; ...
        0, 6*E*I/L^2, 2*E*I/L, 0, -6*E*I/L^2, 4*E*I/L ...
    ];
    T = [ ...
        c, s, 0, 0, 0, 0; ...
        -s, c, 0, 0, 0, 0; ...
        0, 0, 1, 0, 0, 0; ...
        0, 0, 0, c, s, 0; ...
        0, 0, 0, -s, c, 0; ...
        0, 0, 0, 0, 0, 1 ...
    ];
    ke = T' * kl * T;
end

function id = node_id(level, col)
    id = level * 3 + col + 1;
end

function dofs = node_dofs(nodeId)
    base = (nodeId - 1) * 3;
    dofs = [base + 1, base + 2, base + 3];
end

function d = dof_id(level, col, comp)
    d = (node_id(level, col) - 1) * 3 + comp;
end

function write_history_csv(path, time, ag, roof, story, isoMean, baseShear)
    fid = fopen(path, 'w');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, 'time_s,ground_accel_m_per_s2,roof_displacement_m,story_1_drift_ratio,story_2_drift_ratio,story_3_drift_ratio,isolation_displacement_m,total_base_shear_N\n');
    for i = 1:numel(time)
        fprintf(fid, '%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n', ...
            time(i), ag(i), roof(i), story(i, 1), story(i, 2), story(i, 3), isoMean(i), baseShear(i));
    end
end

function write_hysteresis_csv(path, time, isoDisp, isoForce, isoZ)
    fid = fopen(path, 'w');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, 'time_s,iso1_disp_m,iso1_force_N,iso1_z_m,iso2_disp_m,iso2_force_N,iso2_z_m,iso3_disp_m,iso3_force_N,iso3_z_m\n');
    for i = 1:numel(time)
        fprintf(fid, '%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n', ...
            time(i), isoDisp(i, 1), isoForce(i, 1), isoZ(i, 1), isoDisp(i, 2), isoForce(i, 2), isoZ(i, 2), isoDisp(i, 3), isoForce(i, 3), isoZ(i, 3));
    end
end

function write_vector(path, values)
    fid = fopen(path, 'w');
    cleanup = onCleanup(@() fclose(fid));
    for i = 1:numel(values)
        fprintf(fid, '%.12e\n', values(i));
    end
end

function write_json(path, payload)
    try
        text = jsonencode(payload, 'PrettyPrint', true);
    catch
        text = jsonencode(payload);
    end
    fid = fopen(path, 'w');
    cleanup = onCleanup(@() fclose(fid));
    fwrite(fid, text, 'char');
end

function payload = unavailable_payload(status, message)
    payload.schema_version = 'example6-response-v1';
    payload.software = 'MATLAB';
    payload.status = status;
    payload.generated_at_utc = utc_now_string();
    payload.message = message;
    payload.assumptions = { ...
        'No engineering constants were invented; this diagnostic output is written only because the runtime or required input is unavailable.' ...
    };
end

function s = utc_now_string()
    try
        s = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
    catch
        s = datestr(now, 'yyyy-mm-ddTHH:MM:SSZ');
    end
end
