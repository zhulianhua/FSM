% test 1D fft transformation 
N = 32;
x = linspace(-6,6,N);
f = exp(-x.*x/2)/sqrt(2*pi);
y = fft(fftshift(f));
ff = fftshift(ifft(y));
plot(x,f);
hold on;
plot(x,ff,'ro');