function [Rs, stiffness] = computeRsStiffness(lastRs, lastDisp, disp, fy, k0, b)
    % Number of stories
    n = length(lastRs);
    
    % Initialize vectors
    lastShear = zeros(n, 1);
    lastRelativeDisp = zeros(n, 1);
    relativeDisp = zeros(n, 1);
    shear = zeros(n, 1);
    stiffness = zeros(n, 1);
    
    % Calculate shear distribution for the previous step
    for i = 1:n
        lastShear(i) = sum(lastRs(i:n));
    end
    
    % Calculate relative displacements
    lastRelativeDisp(1) = lastDisp(1);
    relativeDisp(1) = disp(1);
    for i = 2:n
        lastRelativeDisp(i) = lastDisp(i) - lastDisp(i-1);
        relativeDisp(i) = disp(i) - disp(i-1);
    end
    
    % Calculate tangent stiffness and shear for each story
    for i = 1:n
        fy_minus_b = fy(i) * (1 - b(i));
        k_sh = b(i) * k0(i);
        deltaDisp = relativeDisp(i) - lastRelativeDisp(i);
        
        c = lastShear(i) + k0(i) * deltaDisp;
        shear(i) = max(k_sh * relativeDisp(i) - fy_minus_b, ...
                       min(k_sh * relativeDisp(i) + fy_minus_b, c));
        
        if abs(shear(i) - c) < 1e-3
            stiffness(i) = k0(i);
        else
            stiffness(i) = k_sh;
        end
    end
    
    % Calculate restoring force
    Rs = zeros(n, 1);
    Rs(n) = shear(n);
    for i = 1:n-1
        Rs(i) = shear(i) - shear(i+1);
    end
end
