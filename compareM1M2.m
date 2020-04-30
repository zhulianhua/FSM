% Compare mouhotFastAlgorithmsComputing2006 with gambaFastSpectralMethod2017
clc; clf;
plot(vi, QFSMHS(:,N/2,N/2),'-x', 'linewidth', 1); % gambaFastSpectralMethod2017
hold on;
plot(vi, QM4(:,N/2,N/2), 's', 'linewidth', 1); % gambaFastSpectralMethod2017
plot(vi, QM10(:,N/2,N/2),'d', 'linewidth', 1); % gambaFastSpectralMethod2017
plot(vi, QM16(:,N/2,N/2),'o', 'linewidth', 1); % gambaFastSpectralMethod2017
xlabel('v')
ylabel('Q(:,N/2,N/2)')
legend({'gamba17, M=14', 'mouhot06 M1=M2=4','mouhot06 M1=M2=10','mouhot06 M1=M2=16'},'FontSize',14)
title('Compare mouhot06 and gamba17, HS model, N = 32')