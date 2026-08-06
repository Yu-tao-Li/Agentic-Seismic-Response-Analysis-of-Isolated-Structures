function K = computeStiffnessMatrix(storyStiffness)
% Assemble the shear-building tangent stiffness matrix.

n = length(storyStiffness);
K = zeros(n, n);
for i = 1:n
    direction = zeros(n, 1);
    direction(i) = 1;
    if i > 1
        direction(i-1) = -1;
    end
    K = K + storyStiffness(i) * (direction * direction');
end
end
