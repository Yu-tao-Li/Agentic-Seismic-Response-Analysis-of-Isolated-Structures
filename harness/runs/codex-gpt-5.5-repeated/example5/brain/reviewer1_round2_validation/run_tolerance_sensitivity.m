function run_tolerance_sensitivity()
% Sweep all three dimensionless tolerances around the adopted values.

scriptDir = fileparts(mfilename('fullpath'));
forceTolerances = [1.0e-6, 1.0e-7, 1.0e-8];
displacementTolerances = [1.0e-8, 1.0e-10, 1.0e-12];
stateTolerances = [1.0e-8, 1.0e-10, 1.0e-12];

adoptedForce = 1.0e-7;
adoptedDisplacement = 1.0e-10;
adoptedState = 1.0e-10;
baseline = Example5_3(adoptedForce, adoptedDisplacement, adoptedState, true);

nRuns = numel(forceTolerances) * numel(displacementTolerances) * numel(stateTolerances);
forceRtol = zeros(nRuns, 1);
displacementRtol = zeros(nRuns, 1);
stateRtol = zeros(nRuns, 1);
peakTopDisplacement = zeros(nRuns, 1);
peakTopAcceleration = zeros(nRuns, 1);
peakIsolationDisplacement = zeros(nRuns, 1);
peakBaseShear = zeros(nRuns, 1);
hysteresisArea = zeros(nRuns, 1);
topDisplacementNrmse = zeros(nRuns, 1);
topAccelerationNrmse = zeros(nRuns, 1);
isolationDisplacementNrmse = zeros(nRuns, 1);
baseShearNrmse = zeros(nRuns, 1);
peakTopDisplacementRelativeChange = zeros(nRuns, 1);
peakTopAccelerationRelativeChange = zeros(nRuns, 1);
peakIsolationDisplacementRelativeChange = zeros(nRuns, 1);
peakBaseShearRelativeChange = zeros(nRuns, 1);
hysteresisAreaRelativeChange = zeros(nRuns, 1);
failedSteps = zeros(nRuns, 1);
maxIterations = zeros(nRuns, 1);
maxFinalResidualRatio = zeros(nRuns, 1);
maxFinalCorrectionRatio = zeros(nRuns, 1);

row = 0;
for i = 1:numel(forceTolerances)
    for j = 1:numel(displacementTolerances)
        for k = 1:numel(stateTolerances)
            row = row + 1;
            current = Example5_3(forceTolerances(i), displacementTolerances(j), ...
                stateTolerances(k), false);

            forceRtol(row) = forceTolerances(i);
            displacementRtol(row) = displacementTolerances(j);
            stateRtol(row) = stateTolerances(k);
            peakTopDisplacement(row) = current.peaks.top_displacement_abs_max;
            peakTopAcceleration(row) = current.peaks.top_acceleration_abs_max;
            peakIsolationDisplacement(row) = current.peaks.isolation_displacement_abs_max;
            peakBaseShear(row) = current.peaks.base_shear_abs_max;
            hysteresisArea(row) = current.peaks.hysteresis_loop_area_abs;

            topDisplacementNrmse(row) = curve_nrmse(current.series.top_displacement, ...
                baseline.series.top_displacement);
            topAccelerationNrmse(row) = curve_nrmse(current.series.top_acceleration, ...
                baseline.series.top_acceleration);
            isolationDisplacementNrmse(row) = curve_nrmse(current.series.isolation_displacement, ...
                baseline.series.isolation_displacement);
            baseShearNrmse(row) = curve_nrmse(current.series.base_shear, ...
                baseline.series.base_shear);

            peakTopDisplacementRelativeChange(row) = relative_change( ...
                current.peaks.top_displacement_abs_max, baseline.peaks.top_displacement_abs_max);
            peakTopAccelerationRelativeChange(row) = relative_change( ...
                current.peaks.top_acceleration_abs_max, baseline.peaks.top_acceleration_abs_max);
            peakIsolationDisplacementRelativeChange(row) = relative_change( ...
                current.peaks.isolation_displacement_abs_max, ...
                baseline.peaks.isolation_displacement_abs_max);
            peakBaseShearRelativeChange(row) = relative_change( ...
                current.peaks.base_shear_abs_max, baseline.peaks.base_shear_abs_max);
            hysteresisAreaRelativeChange(row) = relative_change( ...
                current.peaks.hysteresis_loop_area_abs, baseline.peaks.hysteresis_loop_area_abs);
            failedSteps(row) = current.convergence.failed_step_count;
            maxIterations(row) = current.convergence.max_iterations;
            maxFinalResidualRatio(row) = current.convergence.max_final_residual_ratio;
            maxFinalCorrectionRatio(row) = current.convergence.max_final_correction_ratio;
        end
    end
end

results = table(forceRtol, displacementRtol, stateRtol, ...
    peakTopDisplacement, peakTopAcceleration, peakIsolationDisplacement, ...
    peakBaseShear, hysteresisArea, topDisplacementNrmse, topAccelerationNrmse, ...
    isolationDisplacementNrmse, baseShearNrmse, ...
    peakTopDisplacementRelativeChange, peakTopAccelerationRelativeChange, ...
    peakIsolationDisplacementRelativeChange, peakBaseShearRelativeChange, ...
    hysteresisAreaRelativeChange, failedSteps, maxIterations, ...
    maxFinalResidualRatio, maxFinalCorrectionRatio);
writetable(results, fullfile(scriptDir, 'tolerance_sensitivity_results.csv'));

summary = struct();
summary.run_count = nRuns;
summary.adopted = struct('force_rtol', adoptedForce, ...
    'displacement_rtol', adoptedDisplacement, 'state_rtol', adoptedState);
summary.sweep = struct('force_rtol', forceTolerances, ...
    'displacement_rtol', displacementTolerances, 'state_rtol', stateTolerances);
summary.maximum_over_sweep = struct( ...
    'top_displacement_peak_relative_change', max(peakTopDisplacementRelativeChange), ...
    'top_acceleration_peak_relative_change', max(peakTopAccelerationRelativeChange), ...
    'isolation_displacement_peak_relative_change', max(peakIsolationDisplacementRelativeChange), ...
    'base_shear_peak_relative_change', max(peakBaseShearRelativeChange), ...
    'hysteresis_area_relative_change', max(hysteresisAreaRelativeChange), ...
    'top_displacement_nrmse', max(topDisplacementNrmse), ...
    'top_acceleration_nrmse', max(topAccelerationNrmse), ...
    'isolation_displacement_nrmse', max(isolationDisplacementNrmse), ...
    'base_shear_nrmse', max(baseShearNrmse), ...
    'failed_steps', max(failedSteps), 'max_iterations', max(maxIterations));
summary.baseline_peaks = baseline.peaks;

fid = fopen(fullfile(scriptDir, 'tolerance_sensitivity_summary.json'), 'w');
if fid < 0
    error('Could not write tolerance_sensitivity_summary.json.');
end
fwrite(fid, jsonencode(summary, PrettyPrint=true), 'char');
fclose(fid);

fid = fopen(fullfile(scriptDir, 'tolerance_sensitivity_report.md'), 'w');
if fid < 0
    error('Could not write tolerance_sensitivity_report.md.');
end
fprintf(fid, '# Dimensionless Tolerance Sensitivity\n\n');
fprintf(fid, 'The adopted force-residual, displacement-correction, and state-transition tolerances are `1e-7`, `1e-10`, and `1e-10`, respectively. ');
fprintf(fid, 'A full factorial sweep of 27 combinations varied each criterion over the values listed below.\n\n');
fprintf(fid, '- Force residual tolerance: `1e-6`, `1e-7`, `1e-8`\n');
fprintf(fid, '- Displacement correction tolerance: `1e-8`, `1e-10`, `1e-12`\n');
fprintf(fid, '- State-transition tolerance: `1e-8`, `1e-10`, `1e-12`\n\n');
fprintf(fid, 'Maximum changes relative to the adopted case across all 27 combinations:\n\n');
fprintf(fid, '| Quantity | Maximum relative change or NRMSE |\n');
fprintf(fid, '| --- | ---: |\n');
fprintf(fid, '| Peak top displacement | %.6g |\n', summary.maximum_over_sweep.top_displacement_peak_relative_change);
fprintf(fid, '| Peak top acceleration | %.6g |\n', summary.maximum_over_sweep.top_acceleration_peak_relative_change);
fprintf(fid, '| Peak isolation displacement | %.6g |\n', summary.maximum_over_sweep.isolation_displacement_peak_relative_change);
fprintf(fid, '| Peak base shear | %.6g |\n', summary.maximum_over_sweep.base_shear_peak_relative_change);
fprintf(fid, '| Hysteresis loop area | %.6g |\n', summary.maximum_over_sweep.hysteresis_area_relative_change);
fprintf(fid, '| Top displacement curve NRMSE | %.6g |\n', summary.maximum_over_sweep.top_displacement_nrmse);
fprintf(fid, '| Top acceleration curve NRMSE | %.6g |\n', summary.maximum_over_sweep.top_acceleration_nrmse);
fprintf(fid, '| Isolation displacement curve NRMSE | %.6g |\n', summary.maximum_over_sweep.isolation_displacement_nrmse);
fprintf(fid, '| Base shear curve NRMSE | %.6g |\n\n', summary.maximum_over_sweep.base_shear_nrmse);
fprintf(fid, 'All runs had at most %d failed steps and required at most %d Newton iterations.\n', ...
    summary.maximum_over_sweep.failed_steps, summary.maximum_over_sweep.max_iterations);
fclose(fid);

fprintf('Tolerance sensitivity complete: %d runs.\n', nRuns);
end

function value = relative_change(current, reference)
value = abs(current - reference) / max(abs(reference), eps);
end

function value = curve_nrmse(current, reference)
scale = max(reference) - min(reference);
value = sqrt(mean((current - reference).^2)) / max(scale, eps);
end
