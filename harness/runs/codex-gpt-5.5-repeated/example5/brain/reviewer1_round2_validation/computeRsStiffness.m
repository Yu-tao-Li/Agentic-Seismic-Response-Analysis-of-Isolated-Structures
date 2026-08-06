function [Rs, stiffness, stateRatio] = computeRsStiffness( ...
    lastRs, lastDisp, disp, fy, k0, b, stateRtol)
% Bilinear story-force update using a dimensionless state criterion.

n = length(lastRs);
lastShear = zeros(n, 1);
lastRelativeDisp = zeros(n, 1);
relativeDisp = zeros(n, 1);
shear = zeros(n, 1);
stiffness = zeros(n, 1);
stateRatio = zeros(n, 1);

for i = 1:n
    lastShear(i) = sum(lastRs(i:n));
end

lastRelativeDisp(1) = lastDisp(1);
relativeDisp(1) = disp(1);
for i = 2:n
    lastRelativeDisp(i) = lastDisp(i) - lastDisp(i-1);
    relativeDisp(i) = disp(i) - disp(i-1);
end

for i = 1:n
    fyMinusB = fy(i) * (1 - b(i));
    postYieldStiffness = b(i) * k0(i);
    deltaDisp = relativeDisp(i) - lastRelativeDisp(i);
    elasticTrial = lastShear(i) + k0(i) * deltaDisp;
    shear(i) = max(postYieldStiffness * relativeDisp(i) - fyMinusB, ...
        min(postYieldStiffness * relativeDisp(i) + fyMinusB, elasticTrial));

    forceScale = max([abs(shear(i)), abs(elasticTrial), fy(i), eps]);
    stateRatio(i) = abs(shear(i) - elasticTrial) / forceScale;
    if stateRatio(i) <= stateRtol
        stiffness(i) = k0(i);
    else
        stiffness(i) = postYieldStiffness;
    end
end

Rs = zeros(n, 1);
Rs(n) = shear(n);
for i = 1:n-1
    Rs(i) = shear(i) - shear(i+1);
end
end
