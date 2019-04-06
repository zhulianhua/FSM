% The Fast spectral algorithm from gambaFastSpectralMethod2017
%   
% Reproduce the results of BKW solution for
% for Maxwell-type interactions (gamma = 0)
% NOTE: Only for Maxwell molecules, we can derive an exact 
% solution BKW for the spatial homogeneous Boltzmann equation

% Date: Apr 1, 2019
% Author: Lianhua Zhu
% Depends on lgwt.m and getLebedevSphere.m

clc;clearvars;clf;
tic
%% Computational paramters
% Velocity grid size and number of node for the Guass-Lebelev quadrature
N = 32; M = 14;
% Radial extend for radial integration and velocity domain
R = 6.0; L = (3+sqrt(2.0))*R/4; dv = 2*L/N;
% VHS constant in eq. 1.5, and Page 667
bgamma = 1/4/pi; gamma = 0.0; % 0.0 coressponds to Maxwell molecular model

%% Constant aux arrays for vectorization
li = -N/2:N/2-1;
ll = zeros(3,N,N,N);
vSqr = zeros(N,N,N);
lNorm = zeros(N,N,N);
vi = linspace(-L,L-dv,N);
for l1=1:N
    for l2=1:N
        for l3=1:N
            ll(:,l1,l2,l3) = [li(l1),li(l2),li(l3)]';
            vSqr(l1,l2,l3) = vi(l1)*vi(l1) + vi(l2)*vi(l2) + vi(l3)*vi(l3);
            lNorm(l1,l2,l3) = norm([l1-N/2-1,l2-N/2-1,l3-N/2-1]);
        end
    end
end

%% Distribution function

% allocate data for the velocity distribution function (f)
f = zeros(N,N,N);
% and analytical solution of the Q
Qana = zeros(N,N,N);

% initialize the distribution function (t = 6.5)
% compare the Q(f) with its analytical solution only
t = 6.5;
% Consider to vectorize the loops below
K = 1 - exp(-t/6); Kp = exp(-t/6)/6;
f(:) = 1/(2*(2*pi*K)^1.5)*exp(-vSqr(:)/2/K).*((5*K-3)/K + (1-K)/K/K*vSqr(:));
% analytical solution of Q(f)
Qana(:) = ( (-1.5/K + vSqr(:)/2/K/K).*f(:) + 1/(2*(2*pi*K)^1.5)...
            *exp(-vSqr(:)/2/K).*(3/K/K + (K-2)/K^3*vSqr(:)))*Kp;  

%% Gauss-Legendre and Gauss-Lebedev quadrature points and weights
% N-point Gauss-Legendre quadrature
[GLx,GLw] = lgwt(N,-1,1);
% M-point Gauss-Lebodeve quadrature
leb = getLebedevSphere(M);

%% Storage of G(m,m)
Gmm = zeros(N,N,N);
% compute G(m,m) for VHS (eq. 2.6), vectorized
r = 0.5*R + 0.5*R*GLx;
Gmm(:) = 16*pi*pi*bgamma*0.5*R*(GLw'*(r.^(gamma+2).*(sinc(1.0/L*r.*lNorm(:)'))));

% %Unvectorized
% for l1 = 1:N
%     for l2=1:N
%         for l3=1:N
%             r = 0.5*R + 0.5*R*GLx;
%             y = r.^(gamma+2).*(sinc(pi/L*r*norm([l1-N/2-1,l2-N/2-1,l3-N/2-1])));
%             sy = y'*GLw*R*0.5;
%             Gmm(l1,l2,l3) = 16*pi*pi*bgamma*sy;
%         end
%     end
% end

%% Storage for F(r,k)
Fkr = zeros(N,N,N,N);
% Pre compute F(k,r), eq. 3.8 for VHS model, vectorized
% , actually no need for precomputation when using VHS model
% Frk(:,:) =  4*pi*bgamma*r.^(gamma+2).*sinc(pi/L*r.*lNorm(:)'/2.0);

% %Unvectorized, no need to pre-store for VHS model
% for k1=1:N
%     for k2=1:N
%         for k3=1:N
%             for ri=1:N
%                 r = 0.5*R + 0.5*R*GLx(ri);
%                 Fkr(k1,k2,k3,ri) = 4*pi*bgamma*r^(gamma+2)...
%                     *sinc(1.0/L*r*norm([k1-N/2-1,k2-N/2-1,k3-N/2-1])/2.0);
%             end
%         end
%     end
% end

%% constant aux arrays and temperory space for gain term FFT etc.
lw = zeros(N,N,N);
A = zeros(N,N,N);
B = zeros(N,N,N);

%% a single time calculation evolution
toc % begin computation
fF = fftshift(fftn(fftshift(f)));
% Fourier components of the loss term, N^3log(N) operations
QmF = convnfft(Gmm.*fF, fF);

% Fourier components of the gain term, MN^4log(N) operations
QpF = zeros(N,N,N);
for m=1:M % spherical surface quadratures
    omega = [leb.x(m), leb.y(m), leb.z(m)];
    % compute l*omega, vectorized
    lw(:) = omega*ll(:,:);
    for ri=1:N % Radial quadratures
        % scales
        r = 0.5*R + 0.5*R*GLx(ri);
        wr = GLw(ri);
        A = exp( 1i*pi/L*r*lw/2.0).*fF;
        B = exp(-1i*pi/L*r*lw/2.0).*fF;
        % without storage
        Fkri = 4*pi*bgamma*r^(gamma+2)*sinc(1.0/L*r*lNorm/2.0);
        QpF =  QpF + 0.5*R*wr*leb.w(m).*Fkri.*convnfft(A,B);
    end
end

% iFFT of QF, NOTE of N^3, this is due to the difference 
%  between the presentation in the equations and the matlab
%   real representations
Q = fftshift((ifftn(fftshift(QpF - QmF))))/N^3;
% we need take back real Q
Q = real(Q);

toc % end of computation 

% check if correct
plot(vi, Q(:,16,16),'o','linewidth', 1); % analytical solution of Q, i.e., \partial_t f
hold on;
plot(vi, Qana(:,16,16),'-x', 'linewidth', 1); % numerical solution 
% plot(vi, f(:,16,16),'-','linewidth', 1); % 
xlabel('v')
legend('Q_{num}(:,N/2,N/2)', 'Q_{ana}(:,N/2,N/2)')
title('BKW solution of Maxwell molecular, M=6,N=32')
% should be 3.80e-8 for N=32,M=6, 
fprintf('max difference = %e\n', max(max(max(abs(Q - Qana)))));