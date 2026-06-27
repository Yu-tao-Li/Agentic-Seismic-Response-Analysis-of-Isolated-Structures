function K = computeStiffnessMatrix(storyStiffness)
    % Number of stories
    n = length(storyStiffness);
    
    % Initialize stiffness matrix
    K = zeros(n, n);
    
    % Fill the stiffness matrix
    for i = 1:n
        if i < n
            % Diagonal elements
            K(i, i) = storyStiffness(i) + storyStiffness(i + 1);
            % Off-diagonal elements
            K(i, i + 1) = -storyStiffness(i + 1);
            K(i + 1, i) = -storyStiffness(i + 1);
        else
            % Special treatment for the last row and column
            K(i, i) = storyStiffness(i);
        end
    end
end
