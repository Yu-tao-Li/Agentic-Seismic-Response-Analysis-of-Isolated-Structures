clear; clc; close all;

m = 1000;          % kg per story
k = 500e3;         % N/m per story

M = m * eye(3);

K = k * [ 2  -1   0;
         -1   2  -1;
          0  -1   1];

[Phi, Lambda] = eig(K, M);

lambda = diag(Lambda);
[lambda, order] = sort(lambda, 'ascend');
Phi = Phi(:, order);

omega = sqrt(lambda);
periods = 2 * pi ./ omega;

mode_shapes = zeros(3, 3);
for i = 1:3
    phi = Phi(:, i);
    phi = phi / max(abs(phi));

    if phi(end) < 0
        phi = -phi;
    end

    mode_shapes(i, :) = phi.';
end

omega_rounded = round(omega.' * 1000) / 1000;
periods_rounded = round(periods.' * 1000) / 1000;
mode_shapes_rounded = round(mode_shapes * 1000) / 1000;

results = struct();
results.circular_frequencies_rad_per_s = omega.';
results.periods_s = periods.';
results.normalized_mode_shapes = mode_shapes;
results.circular_frequencies_rad_per_s_rounded = omega_rounded;
results.periods_s_rounded = periods_rounded;
results.normalized_mode_shapes_rounded = mode_shapes_rounded;

fid = fopen('modal_results.json', 'w');
if fid == -1
    error('Could not open modal_results.json for writing.');
end
fprintf(fid, '%s', jsonencode(results));
fclose(fid);

disp('Circular frequencies (rad/s):');
disp(omega_rounded);

disp('Periods (s):');
disp(periods_rounded);

disp('Normalized mode shapes, rows = modes, columns = stories 1 to 3:');
disp(mode_shapes_rounded);

stories = [0 1 2 3];

for i = 1:3
    amplitudes = [0 mode_shapes_rounded(i, :)];

    figure('Color', 'w');
    plot(amplitudes, stories, '-o', ...
        'LineWidth', 2, ...
        'MarkerSize', 7, ...
        'MarkerFaceColor', [0.1 0.35 0.75]);
    grid on;
    xlabel('Normalized amplitude');
    ylabel('Story number');
    title(sprintf('Mode %d Shape', i));
    xlim([-1.1 1.1]);
    ylim([0 3]);
    yticks(0:3);
    yticklabels({'Ground', '1', '2', '3'});
    xline(0, '--k');
    saveas(gcf, sprintf('mode_%d_shape.png', i));
end
