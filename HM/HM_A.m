%HM FFT
 sigma_0 = 1/32;
 Q = 1;
 L = 1;
 N = 256;
 dx = L/N;
 gauss     = @(x) Q/sqrt(pi*sigma_0^2)*exp(-( ((x-N/2*dx).^2) / sigma_0^2));
 gauss_hat = @(k) Q*exp(-( ((k).^2) * sigma_0^2)/4);
 X1 = linspace(1,N,N);
 A1 = gauss(X1*dx);
 B1 = fft(A1)/N;
  k = linspace(0,N/2,N/2);
 B2 = gauss_hat(k); 
 loglog(k,abs(B1(1:N/2)),k,B2);
 title('Fourier Coefficient N=256');
 legend('Coefficient from the FFT','Coefficient from the unbound domain')

%Y_hat = gauss_hat(1);
%plot(K,abs(Y_fft(1:length(K))/N),K,Y_hat);
% Y_ex = gauss(X);
% Y_re = zeros(1,N);
% for m = 1 : N
%     for n = 1 : length(Y_fft)/2
%        Y_re(m) = Y_re(m) + Y_fft(n) * exp(i*k(n)*X(m))/N;
%     end
%     Y_re(m) = real(F_hat(1))/N;
%     for k = 1 : length(F_hat)/2-1
%         Y_re(m) = 2*( real(F_hat(k))*cos(k*X(m)) - imag(F_hat(k))*sin(k*X(m)) )/N + Y_re(m);  
%     end
%     Y_re(m) = 2*( real(F_hat(N/2))*cos((N/2)*X(m) ) )/N + Y_re(m);
% end
%plot(X,Y_re),X,Y_ex);    
% Y_HAT = 0;
% for n = 1 : length(Y_fft)
%     Y_HAT = Y_fft(n)*exp(-i*K(n)*
% end

%Analysis
% k_starhF = @(kh) (exp(-i*2*kh)-6*exp(-i*kh)+3+2*exp(i*kh))/6;
% kh = linspace(0,pi,100);
% kh_star = abs(k_starhF(kh));
% khpi = kh/pi;
% K_A1 = khpi;
% K_A2 = kh_star/pi;
% K_A3 = sin(kh)/pi;
% K_A4 = (4/3*sin(kh)-1/6*sin(2*kh))/pi;
% K_A5 = (45/30*sin(kh) - 9/30*sin(2*kh) + 1/30*sin(3*kh))/pi;
% L1 = ones(1,100)*0.25;
% L2 = ones(1,100)*0.5;
% L3 = ones(1,100)*0.75;
% figure;
% plot(khpi,K_A1,khpi,K_A2,khpi,K_A3,khpi,K_A4,khpi,K_A5,L1,khpi,'Black',L2,khpi,'Black',L3,khpi,'Black');
% legend('Ideal Case','Order 4','E2','E4','E6');
% figure;
% plot(khpi,K_A1./khpi,khpi,K_A2./khpi,khpi,K_A3./khpi,khpi,K_A4./khpi,khpi,K_A5./khpi);
% legend('Ideal Case','Order 4','E2','E4','E6');

% Lambda
% Stability E2




%h_lambda_E2 = @(kh) -( 1-cos(kh)+i*sin(kh));
%h_lambda_E2 = @(kh) -( 1-exp(-i*kh));
%h_lambda_E2 = @(kh) -( 1-exp(-i*kh));
% h_lambda_E4 = @(kh) -( 4/6 * (exp(i*kh)-exp(-i*kh)) - 1/12*(exp(i*2*kh)-exp(-i*2*kh)) );
% 
% k_starhF = @(kh) (exp(-i*2*kh)-6*exp(-i*kh)+3+2*exp(i*kh))/6;
% h_lambda_S3 = @(kh) -(exp(-i*2*kh)-6*exp(-i*kh)+3+2*exp(i*kh))/6;
% k_star_h_E2 = @(kh) -i*(sin(kh));
% k_star_h_E4 = @(kh) -i*(8/6*sin(kh) - 1/6 * sin(2*kh));
% k_star_h_E6 = @(kh) -i*(45/30*sin(kh) - 9/30 * sin(2*kh) + 1/30*sin(3*kh));

% E2X = linspace(-4,4,100);
% E2Y = h_lambda_E4(E2X); 
% 
% 
% S3X = linspace(-4,4,100);
% S3Y = h_lambda_S3(S3X);
% plot(real(S3Y),imag(S3Y),real(E2Y),imag(E2Y));
% axis('equal');

% figure;
% Q3X = linspace(0,pi,100);
% Q3Y = h_lambda_S3(Q3X);
% plot(Q3X/pi,-angle(h_lambda_S3(Q3X))/pi,'--',Q3X/pi,-angle(k_star_h_E2(Q3X))/pi,Q3X/pi,-angle(k_star_h_E4(Q3X))/pi,'.',Q3X/pi,-angle(k_star_h_E6(Q3X))/pi,Q3X/pi,angle(Q3X*1/pi),'black--');
% xlabel('kh{\pi}');
% ylabel('h\lambda_g/ c \pi');
% title('Phase Error');
% legend('Decentered O3','Centered E2','Centered E4','Centered E6','Exact Solution');


% Q3X = linspace(0,pi,50);
% Q3S3 = h_lambda_S3(Q3X);
% Q3E2 = k_star_h_E2(Q3X);
% Q3E4 = k_star_h_E4(Q3X);
% Q3E6 = k_star_h_E6(Q3X);
% plot(real(Q3S3),imag(Q3S3),real(Q3E2),imag(Q3E2),'o',real(Q3E4),imag(Q3E4),'+',real(Q3E6),imag(Q3E6),'*');
% title('Stability Curve in complex plane');
% xlabel('Real \lambda_k');
% ylabel('Imag \lambda_k');
% legend('Decentered O3','Centered E2','Centered E4','Centered E6');



% RK4XV = linspace(-4,3,100);
% RK4YV = linspace(-4,3,100);
% [RK4X,RK4Y] = meshgrid(RK4XV,RK4YV);
% RK4Z = RK4X + i* RK4Y ; 
% RK4F = 1 + RK4Z + 1/2 *(RK4Z).^2 + 1/6 * RK4Z.^3 + 1/24 * RK4Z.^4 ;
% RK4FMAG = abs(RK4F);
% hold on;
% contour(RK4X,RK4Y,RK4FMAG,[1 1]);
% axis('equal')
% hold off;



% AFUN = @(kh) sin(kh);
% A = linspace(-4,1,100);
% B = linspace(-4,1,100);
% X = A + i*B;
% Y = AFUN(X);
% plot(real(Y),imag(Y));
% axis('equal'); 
