% The Direct spectral algorithm from gambaFastSpectralMethod2017
%   
% Attempting to reproduce the results of BKW solution for
% for Maxwell-type interactions (gamma = 0)
% NOTE: Only for Maxwell molecules, we can derive an exact 
% solution BKW for the spatial homogeneous Boltzmann equation

% Date: Mar 24, 2019
% Author: Lianhua Zhu

% Velocity grid size
N = 32;

% Radial extend for radial integration
R = 6.0;

% Velocity domain
L = (3+sqrt(2.0))*R/4;
% Uniform velocity grid spacing
dv = 2*L/N;
% Allocate the UNIFORM velocity grid in [-L, L]^3
vi = linspace(-L+dv,L,N);
[vx, vy, vz] = meshgrid(vi,vi,vi);

% VHS constant in eq. 1.5
bgamma = 0.0;
gamma = 0.0; % 0.0 coressponds to Maxwell molecular model

% allocate data for the velocity distribution function (f) and
% its Fourier space components (fF)
f = zeros(N,N,N);
% initialize the distribution function (t = 6.5)
% compare the Q(f) with its analytical solution only
t = 6.5;
% Consider to vectorize the loops below
for l1 = 1:N
    for l2 = 1:N
        for l3 = 1:N
            K = 1 - exp(-t/6);
            vSqr = vi(l1)*vi(l1) + vi(l2)*vi(l2) + vi(l3)*vi(l3);
            f(l1,l2,l3) = 1/(2*(2*pi*K)^1.5) ...
                *exp(-vSqr/2/K)*((5*K-3)/K + (1-K)/K/K*vSqr);
        end
    end
end

fF = fftshift(ifftn(fftshift(f)));
% Used for checking fftn and ifftn
% ff = fftshift(fftn(fftshift(fF)));

% N-point Gauss-Legendre quadrature
GLxw = load('GLxwN32.dat');
GLx = GLxw(:,1);
GLw = GLxw(:,2);

% % Storage of G(m,m)
% G = zeros(N,N,N); % m is 3-componnet index
% % compute G(m,m) for VHS (eq. 2.6)
% for l1 = 1:N
%     for l2=1:N
%         for l3=1:N
%             y = zeros(N,1);
%             for ri = 1:N % radial direction
%                 % m = k- N/2 -1;
%                 r = GLx(ri);
%                 xx = pi/L*r*norm([l1-N/2-1,l2-N/2-1,l3-N/2-1]); % l = m, so |l+m|/2 = |m|
%                 y(ri) = r^(gamma+2)*(sin(xx)/xx);
%             end
%             sy = y*GLw';
%             G(l1,l2,l3) = 16*pi*pi*bgamma*sy;
%         end
%     end
% end

% Storage of G(l,m), 6 dimensional, storage N^6
Glm = zeros(N,N,N,N,N,N);
for l1=1:N
    fprintf('i = %d\n', l1);
    for l2=1:N
        for l3=1:N
            fprintf('l1 = %d, l2 = %d, l3 = %d\n',l1, l2, l3);                
            for m1=1:N
                for m2=1:N
                    for m3=1:N
                        y = zeros(N,1);
                        for ri = 1:N
                            % scale from [-1,1] to [0, R] of GL 
                            % [mathewsGaussLegendreIntegration]
                            r = 0.5*R + 0.5*R*GLx(ri);
                            ll1 = l1-N/2-1;
                            ll2 = l2-N/2-1;
                            ll3 = l3-N/2-1;
                            mm1 = m1-N/2-1;
                            mm2 = m2-N/2-1;
                            mm3 = m3-N/2-1;
                            xxp = pi/L*r*norm([ll1,ll2,ll3] ...
                                  + [mm1,mm2,mm3])/2.0; % l = m, so |l+m|/2 = |m|
                            xxm = pi/L*r*norm([ll1,ll2,ll3] ...
                                  - [mm1,mm2,mm3])/2.0;
                            y(ri) = r^(gamma+2)*(sin(xxp)/xxp)*(sin(xxm)/xxm);
                        end
                        sy = y'*GLw*R*0.5; % scale back
                        Glm(l1,l2,l3,m1,m2,m3) = 16*pi*pi*bgamma*sy;
                    end
                end
            end
        end
    end
end

save('GlmN32',Glm);

% time evolution
Q = zeros(N,N,N);
QF = zeros(N,N,N);
for t = 1:1
    % FFT
    for k1=1:N
        for k2=1:N
            for k3=1:N
                QF(k1,k2,k3) = 0.0;
                % convolution
                for m1=1:N
                    for m2=1:N
                        for m3=1:N
                            l1 = k1 - m1;
                            l2 = k2 - m2;
                            l3 = k3 - m3;
                            QF(k1,k2,k3) = ...
                                QF(k1,k2,k3) + ...
                                (G(l1,l2,l3,m1,m2,m3) - G(m1,m2,m3,m1,m2,m3))...
                                *fF(l1,l2,l3)*fF(m1,m2,m3);
                        end
                    end
                end
            end
        end
    end
    % iFFT
    Q = fftshift(fftn(fftshift(QF)));
    % update f
    % f = f + Q*dt;
end