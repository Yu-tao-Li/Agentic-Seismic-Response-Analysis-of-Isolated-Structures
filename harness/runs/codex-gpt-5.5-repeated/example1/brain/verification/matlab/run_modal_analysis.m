% Modal analysis of a three-story shear building in SI units.
% Rows in mode shape output are floors 1..3; columns are modes 1..3.

clear; clc;

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir)
    this_dir = pwd;
end

m_story = 1000.0;       % kg
k_story = 500000.0;    % N/m

M = m_story * eye(3);
K = k_story * [ 2.0, -1.0,  0.0; ...
               -1.0,  2.0, -1.0; ...
                0.0, -1.0,  1.0 ];

[V, D] = eig(K, M);
lambda = real(diag(D));
[lambda_sorted, order] = sort(lambda, 'ascend');
mode_shapes = real(V(:, order));
circular_frequencies = sqrt(lambda_sorted);
periods = 2.0 * pi ./ circular_frequencies;

for mode_id = 1:size(mode_shapes, 2)
    mode_shapes(:, mode_id) = normalize_mode(mode_shapes(:, mode_id));
end

results = struct();
results.software = 'MATLAB';
results.execution = struct( ...
    'matlab_version', version, ...
    'source_file', fullfile(this_dir, 'run_modal_analysis.m'));
results.model = struct( ...
    'description', 'Three-story shear building with equal lumped floor masses and equal story stiffnesses.', ...
    'units', 'SI: kg, N, m, s', ...
    'story_mass_kg', m_story, ...
    'story_stiffness_N_per_m', k_story, ...
    'mass_matrix_kg', M, ...
    'stiffness_matrix_N_per_m', K);
results.mode_shape_convention = 'Rows are floors 1..3; columns are modes 1..3; each mode is normalized to max(abs(component)) = 1 and signed so the roof component is positive.';
results.circular_frequencies_rad_per_s = circular_frequencies(:).';
results.periods_s = periods(:).';
results.normalized_mode_shapes = mode_shapes;

json_text = jsonencode(results, "PrettyPrint", true);
json_path = fullfile(this_dir, 'modal_results.json');
fid = fopen(json_path, 'w');
assert(fid > 0, 'Could not open modal_results.json for writing.');
fprintf(fid, '%s\n', json_text);
fclose(fid);

txt_path = fullfile(this_dir, 'modal_results.txt');
fid = fopen(txt_path, 'w');
assert(fid > 0, 'Could not open modal_results.txt for writing.');
fprintf(fid, 'Three-story shear building modal analysis - MATLAB\n');
fprintf(fid, 'Units: kg, N, m, s\n');
fprintf(fid, 'Story mass: %.12g kg\n', m_story);
fprintf(fid, 'Story stiffness: %.12g N/m\n\n', k_story);
fprintf(fid, 'Mass matrix M [kg]:\n');
write_matrix(fid, M);
fprintf(fid, '\nStiffness matrix K [N/m]:\n');
write_matrix(fid, K);
fprintf(fid, '\nMode results:\n');
fprintf(fid, 'Mode, circular_frequency_rad_per_s, period_s, normalized_floor_1, normalized_floor_2, normalized_floor_3\n');
for mode_id = 1:3
    fprintf(fid, '%d, %.15g, %.15g, %.15g, %.15g, %.15g\n', ...
        mode_id, circular_frequencies(mode_id), periods(mode_id), ...
        mode_shapes(1, mode_id), mode_shapes(2, mode_id), mode_shapes(3, mode_id));
end
fclose(fid);

disp('MATLAB modal analysis completed.');
disp(json_path);

function phi_norm = normalize_mode(phi)
    phi_norm = real(phi(:));
    scale = max(abs(phi_norm));
    assert(scale > 0.0, 'Zero mode shape encountered.');
    phi_norm = phi_norm ./ scale;
    if phi_norm(end) < 0.0
        phi_norm = -phi_norm;
    end
end

function write_matrix(fid, A)
    for row_id = 1:size(A, 1)
        fprintf(fid, '  ');
        fprintf(fid, '% .15g ', A(row_id, :));
        fprintf(fid, '\n');
    end
end
