%figure;
clf;
semilogy(ts(2:end), maxErrN8(2:end), 'k-v');hold on;
semilogy(ts(2:end), maxErrN12(2:end),'g-x');
semilogy(ts(2:end), maxErrN16(2:end),'b-s');
semilogy(ts(2:end), maxErrN32(2:end),'r-o');

ylim([1e-8,1e0])
legend('N8', 'N12','N16','N32');
% legend('N12','N16','N32');


title('max relative norm error of f with t');