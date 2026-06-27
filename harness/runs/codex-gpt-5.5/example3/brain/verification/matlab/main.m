clear; clc;
dt = 0.02;
g = 9.81;
targetPGA = 0.40 * g;
raw = load(fullfile('..','..','input_data','GM1.txt'));
raw = raw(:);
scale = targetPGA / max(abs(raw));
ag = raw * scale;
n = length(ag);
t = (0:n-1)' * dt;
masses = [500 1000 1000 1000];
storyK = [200000 500000 500000 500000];
ndof = 4;
M = diag(masses);
K = zeros(ndof);
for spring = 1:ndof
    if spring == 1
        K(1,1) = K(1,1) + storyK(1);
    else
        left = spring - 1;
        right = spring;
        K(left,left) = K(left,left) + storyK(spring);
        K(right,right) = K(right,right) + storyK(spring);
        K(left,right) = K(left,right) - storyK(spring);
        K(right,left) = K(right,left) - storyK(spring);
    end
end
[V,D] = eig(K,M);
omega = sort(sqrt(diag(D)));
xi = 0.05;
alpha = 2*xi*omega(1)*omega(2)/(omega(1)+omega(2));
betaRay = 2*xi/(omega(1)+omega(2));
C = alpha*M + betaRay*K;
C(1,1) = C(1,1) + 2000;
gamma = 0.5;
beta = 0.25;
u = zeros(ndof,n);
v = zeros(ndof,n);
a = zeros(ndof,n);
r = ones(ndof,1);
a(:,1) = M \ (-C*v(:,1) - K*u(:,1) - M*r*ag(1));
Keff = M + gamma*dt*C + beta*dt^2*K;
for i = 1:n-1
    up = u(:,i) + dt*v(:,i) + 0.5*dt^2*(1-2*beta)*a(:,i);
    vp = v(:,i) + dt*(1-gamma)*a(:,i);
    rhs = -M*r*ag(i+1) - C*vp - K*up;
    a(:,i+1) = Keff \ rhs;
    u(:,i+1) = up + beta*dt^2*a(:,i+1);
    v(:,i+1) = vp + gamma*dt*a(:,i+1);
end
topU = u(ndof,:)';
topA = a(ndof,:)';
outDis = [t, topU];
outAcc = [t, topA];
save('topStoDis.txt', 'outDis', '-ascii');
save('topStoAcc.txt', 'outAcc', '-ascii');
resp = struct();
resp.software = 'MATLAB';
resp.time = t;
resp.top_displacement = topU;
resp.top_acceleration = topA;
resp.peak_top_displacement_m = max(abs(topU));
resp.peak_top_acceleration_m_s2 = max(abs(topA));
resp.input_pga_normalized_m_s2 = targetPGA;
resp.raw_pga_m_s2 = max(abs(raw));
resp.scale_factor = scale;
resp.dt = dt;
resp.n_samples = n;
resp.rayleigh_alpha = alpha;
resp.rayleigh_beta = betaRay;
fid = fopen('response_output.json','w');
fprintf(fid, '%s', jsonencode(resp, 'PrettyPrint', true));
fclose(fid);
