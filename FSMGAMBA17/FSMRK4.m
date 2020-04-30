clc;
clf;
% clearvars;

N = 32;
M = 6;
bgamma = 1/4/pi;
gamma = 0;
R = 6.0;

dt = 0.1;
tBegin = 5.5;
tEnd = 8.5;
Nt = (tEnd-tBegin)/dt;
ts = tBegin:dt:tEnd;

f = zeros(N,N,N);
fExt = zeros(N,N,N);

maxErr = zeros(Nt+1,1);

K1 = f;
K2 = K1;
K3 = K1;
K4 = K1;

fsm = FSMClass(6,1/4/pi,0,N,M);

% at the beginning
t = tBegin;
K = 1 - exp(-t/6);
f(:) = 1/(2*(2*pi*K)^1.5)*exp(-fsm.vSqr(:)/2/K).*((5*K-3)/K + (1-K)/K/K*fsm.vSqr(:));
fExt(:) = f(:);

maxErr(1) = 0;
tic
for ti = 1:Nt % from ts(2) to ts(end)
    fprintf("step %d calculating...\n", ti);
    K1 = dt*fsm.getQ(f);
    K2 = dt*fsm.getQ(f+0.5*K1);
    K3 = dt*fsm.getQ(f+0.5*K2);
    K4 = dt*fsm.getQ(f+K3);
    % new time step f
    f = f + (K1 + 2*K2 + 2*K3 + K4)/6.0;  
    % new time step exact f
    t = t + dt;
    % cal analytical solution at the new time step
    K = 1 - exp(-t/6);
    fExt(:) = 1/(2*(2*pi*K)^1.5)*exp(-fsm.vSqr(:)/2/K).*((5*K-3)/K + (1-K)/K/K*fsm.vSqr(:));
    % new time error
    maxErr(ti+1) = norm(f(:)-fExt(:),inf)./norm(fExt(:),inf);
end
toc

%figure;
semilogy(ts(2:end), maxErr(2:end),'r-o');
% plot(f(:,N/2,N/2), 'r-o'); hold on;
% plot(fExt(:,N/2,N/2), 'b-');
ylim([1e-8,1e0])
