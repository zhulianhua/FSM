%% The Fast spectral algorithm of wuDeterministicNumericalSolutions2013
%   
% 3D Maxwell gas model
%  $ B = b_gamma $

% NOTE: Comapare with the analytical BKW solution, use exactly the
%   parameters as in wuDeterministicNumericalSolutions2013
% NOTE: in gambaFastSpectralMethod2017, sqrt{RT} is the velocity scale
%   while in this study, the scale is sqrt{2RT}, the following computation
%   will be non-dimensional form

% Date: Apr 25, 2019
% Author: Lianhua Zhu
% Depends on lgwt.m for Gauss-Legendre quadrature, avaiable online

%% Computational paramters
figure;
clc; hold on;
tic

% Velocity grid size and number of node for the polar and azimuthal angles
N=32;
M1 = 6; M2 = 6;

% Radial extend for radial integration and velocity domain
R = 6.0;
L = 8.0;

dv = 2*L/N;

% (3+sqrt(2.0))*R/2; dv = 2*L/N;
% VHS constant in eq. 1.5, and Page 667
KnPrim = pi*32/5;
% KnPrim = pi;
gamma = 0.0; % 0.0 coressponds to Maxwell model

%% Constant arrays assisting vectorizations later
li = -N/2:N/2-1;
ll = zeros(3,N,N,N); 
vSqr = zeros(N,N,N);
lNorm = zeros(N,N,N);
% how to choose the discrete velocity?
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

%% Initial distribution function
% allocate data for the velocity distribution function (f)
f = zeros(N,N,N);
f1 = zeros(N,N,N);

% initialize the distribution function (t = 6.5), 

% NOTE, below is based on Lei's t
% t0 = 0;
% t1 = 1e-4;
% % Consider to vectorize the loops below
% K0 = 1 - 0.4*exp(-t0/6);
% K1 = 1 - 0.4*exp(-t1/6);
% K = K0;
% f(:) = 1/(2*(2*pi*K)^1.5)*exp(-vSqr(:)/2/K).*((5*K-3)/K + (1-K)/K/K*vSqr(:));
% K = K1;
% f1(:) = 1/(2*(2*pi*K)^1.5)*exp(-vSqr(:)/2/K).*((5*K-3)/K + (1-K)/K/K*vSqr(:));
% Qana = (f1 - f)/(t1-t0);


% NOTE, below is based on gambaFastSpectralMethod2017
t = 6.5;
% Consider to vectorize the loops below

K = 1 - exp(-t/6); Kp = exp(-t/6)/6;
f(:) = 1/(2*(2*pi*K)^1.5)*exp(-vSqr(:)/2/K).*((5*K-3)/K + (1-K)/K/K*vSqr(:));
Qana = zeros(N,N,N);
% analytical solution of Q(f)
Qana(:) = ( (-1.5/K + vSqr(:)/2/K/K).*f(:) + 1/(2*(2*pi*K)^1.5)...
            *exp(-vSqr(:)/2/K).*(3/K/K + (K-2)/K^3*vSqr(:)))*Kp;


%% Gauss-Legendre quadrature points and weights
% N-point Gauss-Legendre quadrature. NOTE: we have used the scaled and shifted x
[GLx,GLw] = lgwt(N,0,R);

%% Storage of G(m,m)
Gmm = zeros(N,N,N);

for p=1:M1
    for q = 1:M2
        thetaPrim = p*pi/M1;
        varphi = q*pi/M2;
        
        e = [sin(thetaPrim)*cos(varphi), ...
             sin(thetaPrim)*sin(varphi), ...
             cos(thetaPrim)...
            ];
        le = 1.0/L * (e * ll(:,:))';
        
        % first we calculate the projection of l on e, (l is m in my note)
        projectionOnE =  vecnorm( (repmat(e,[N^3,1]) .* (e*ll(:,:))')')';
        
        % Then we use the GOUGU theorm to get the norm of the projection 
        %  of the velocity m on the perpenticular plane of the unit vector e
        PiePerpm = pi/L * sqrt(lNorm(:).^2 - projectionOnE.^2);
        % GL quadrature, we have assumed $b(\rho) = 1$
        sumy = zeros(N^3,1);
        for ri = 1:N
            sumy = sumy + GLw(ri) * GLx(ri)*(besselj(0, GLx(ri) * PiePerpm));
        end
        sumy = 2*pi*sumy;
        Gmm(:) = Gmm(:) + 2*R*sinc(R*le).* sumy .* sin(thetaPrim); 
    end
end

%% Storage for alpha'(m,p,q)
alphaPrimpqm = zeros(N,N,N,M1,M2);
% In the following code, vectorized in the first 3 dimension (:,:,:,p,q)
for p=1:M1
    for q = 1:M2
        
        thetaPrim = p*pi/M1;
        varphi = q*pi/M2;
        
        e = [sin(thetaPrim)*cos(varphi), ...
             sin(thetaPrim)*sin(varphi), ...
             cos(thetaPrim)...
            ];
         
        % first we calculate the projection of l on e, (l is m in my note)
        projectionOnE =  vecnorm( (repmat(e,[N^3,1]) .* (e*ll(:,:))')')';
        
        % Then we use the GOUGU theorm to get the norm of the projection 
        %  of the velocity m on the perpenticular plane of the unit vector e
        PiePerpm = pi/L * sqrt(lNorm(:).^2 - projectionOnE.^2);
        
        % GL quadrature, we have assumed $b(\rho) = 1$
        sumy = 0.0;
        for ri = 1:N
            sumy = sumy + GLw(ri) * GLx(ri)*(besselj(0, GLx(ri) * PiePerpm));
        end
        alphaPrimpqm(:,:,:,p,q) = reshape(sumy * 2 * pi, N, N, N);
    end
end

%% Storage for alpha(l,p,q)
%   Can be computed 'on-the-fly', but here we pre-compute and store it
alphapqm = zeros(N,N,N,M1,M2);
for p = 1:M1
    for q = 1:M2
        thetaPrim = p*pi/M1;
        varphi = q*pi/M2;
        
        e = [sin(thetaPrim)*cos(varphi), ...
             sin(thetaPrim)*sin(varphi), ...
             cos(thetaPrim)...
            ];
        
        % Calcualte the projection (norm) of the velocity m on e.
        %   NOTE: the sin(thetaPrim) is mistakenly omitted in 
        %   mouhotFastAlgorithmsComputing2006,
        le = 1.0/L * (e * ll(:,:))';       
        alphapqm(:,:,:,p,q) = reshape(2*R*sinc(R*le),  N, N, N);
    end
end

%% A single time calculation evolution
toc % begin computation
fF = fftshift(fftn(fftshift(f)));
% Fourier components of the loss term, N^3log(N) operations
QmF = convnfft(Gmm.*fF, fF);
QmF = QmF * 4 * pi*pi/M1/M2/KnPrim;
% Fourier components of the gain term, M1*M2*N^3*log(N) operations,
%  opposed to M*N^4*log(N) in gambaFastSpectralMethod2017
QpF = zeros(N,N,N);
for p=1:M1
    for q=1:M2
        thetaPrim = p*pi/M1;
        % scales
        A = alphapqm(:,:,:,p,q).*fF;
        B = alphaPrimpqm(:,:,:,p,q).*fF;
        QpF =  QpF +  sin(thetaPrim)* convnfft(A,B);
    end
end
QpF = QpF * 4* pi*pi/M1/M2/KnPrim;

% iFFT of QF
%   NOTE of the 1/N^3 scale, this is due to the difference between the 
%   presentation in the equations and the matlab's real representations
Q = fftshift((ifftn(fftshift(QpF - QmF))))/N^3;
% we need take back real Q
Q = real(Q);

toc % end of computation

%% check if correct
sum(abs(Q(:)-Qana(:)))/sum(abs(Qana(:)))


plot(vi, Q(:,N/2,N/2),'o','linewidth', 1); % analytical solution of Q, i.e., \partial_t f
hold on;
plot(vi, Qana(:,N/2,N/2),'-x', 'linewidth', 1); % numerical solution 
% plot(vi, f(:,16,16),'-','linewidth', 1); % 
xlabel('v')
legend('Q_{num}(:,N/2,N/2)', 'Q_{ana}(:,N/2,N/2)')
title('BKW solution of Maxwell molecular, M=6,N=32')
% should be 3.80e-8 for N=32,M=6, 
fprintf('max difference = %e\n', max(max(max(abs(Q - Qana)))));