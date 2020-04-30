tic;
N = 10000;
A = zeros(N,N);
B = zeros(N,N);
C = zeros(N,N);

for i=1:N
      for j=1:N
          A = (i*2+j);
          B = (j*2+i);
      end
end

s = 0.0;
for t = 1:100
    At = A*t^2 + 2.0;
    Bt = B*t*2+1;
    C = At./Bt;
    s = s + sum(sum(C));  
end
fprintf("s = %f\n", s);
toc