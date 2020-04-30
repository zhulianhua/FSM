classdef FSMClass
   properties
      % scalars, constractor paramters
      R;
      bgamma;gamma;
      N;M;
      dv;
      % scalars, computed
      L;
      % struct
      leb;
      % M array;
      lebx, leby, lebz, lebw;
      % N array
      GLx; GLw; li; vi;
      % (3,N,N,N) array
      ll;
      % (N,N,N) array
      vSqr; lNorm; Gmm; lw; fF; QmF; QpF;
   end
   methods
      function obj = FSMClass(Rin, bgammain, gammain, Nin, Min)
         if nargin == 5
            obj.R = Rin;
            obj.L = (3+sqrt(2.0))*obj.R/4;
            obj.bgamma = bgammain;
            obj.gamma = gammain;
            obj.N = Nin;
            obj.M = Min;
            
            N = obj.N;
            
            obj.dv = 2*obj.L/obj.N;
            
            obj.li = -obj.N/2:obj.N/2-1;
%             obj.li = -obj.N/2+0.5:obj.N/2-0.5;
            obj.vi = linspace(-obj.L,obj.L-obj.dv,obj.N);
%             obj.vi = linspace(-obj.L+0.5*obj.dv, obj.L-0.5*obj.dv,obj.N);
            
            obj.ll = zeros(3,obj.N,obj.N,obj.N);
            obj.vSqr = zeros(obj.N,obj.N,obj.N);
            obj.lNorm = zeros(obj.N,obj.N,obj.N);
           
            for l1=1:obj.N
                for l2=1:obj.N
                    for l3=1:obj.N
                        obj.ll(:,l1,l2,l3) = [obj.li(l1),obj.li(l2),obj.li(l3)]';
                        obj.vSqr(l1,l2,l3) = obj.vi(l1)*obj.vi(l1) ...
                            + obj.vi(l2)*obj.vi(l2) + obj.vi(l3)*obj.vi(l3);
                        obj.lNorm(l1,l2,l3) = norm([l1-obj.N/2-1,l2-obj.N/2-1,l3-obj.N/2-1]);
%                         obj.lNorm(l1,l2,l3) = norm([l1-obj.N/2-0.5,l2-obj.N/2-0.5,l3-obj.N/2-0.5]);
                    end
                end
            end
            
            [obj.GLx,obj.GLw] = lgwt(obj.N,-1,1);
            % M-point Gauss-Lebodeve quadrature
            leb = getLebedevSphere(obj.M);
            obj.lebx = leb.x;
            obj.leby = leb.y;
            obj.lebz = leb.z;
            obj.lebw = leb.w;
            %% Storage of G(m,m)
            obj.Gmm = zeros(N,N,N);
            obj.lw = zeros(N,N,N);
            obj.fF = zeros(N,N,N);
            obj.QmF = zeros(N,N,N);
            obj.QpF = zeros(N,N,N); 
            
            % compute G(m,m) for VHS (eq. 2.6), vectorized
            r = 0.5*obj.R + 0.5*obj.R*obj.GLx;
            obj.Gmm(:) = 16*pi*pi*obj.bgamma*0.5*obj.R...
                *(obj.GLw'*(r.^(obj.gamma+2).*(sinc(1.0/obj.L*r.*obj.lNorm(:)'))));
         else
             error("Constructor var error");
         end
      end
      
      function Q = getQ(obj,f)
        obj.fF = fftshift(fftn(fftshift(f)));

        % Fourier components of the loss term, N^3log(N) operations
        obj.QmF = convnfft(obj.Gmm.*obj.fF, obj.fF);
        % Fourier components of the gain term, MN^4log(N) operations
        obj.QpF(:) = 0.0;
        for m=1:obj.M % spherical surface quadratures
            omega = [obj.lebx(m), obj.leby(m), obj.lebz(m)];
            % compute l*omega, vectorized
            obj.lw(:) = omega*obj.ll(:,:);
            for ri=1:obj.N % Radial quadratures
                % scales
                r = 0.5*obj.R + 0.5*obj.R*obj.GLx(ri);
                wr = obj.GLw(ri);
                A = exp( 1i*pi/obj.L*r*obj.lw/2.0).*obj.fF;
                B = exp(-1i*pi/obj.L*r*obj.lw/2.0).*obj.fF;
                Fkri = 4*pi*obj.bgamma*r^(obj.gamma+2)*sinc(1.0/obj.L*r*obj.lNorm/2.0);
                obj.QpF = obj.QpF + 0.5*obj.R*wr*obj.lebw(m).*Fkri.*convnfft(A,B);
            end
        end
        Q = fftshift((ifftn(fftshift(obj.QpF - obj.QmF))))/obj.N^3;
        Q = real(Q);
      end
   end
end