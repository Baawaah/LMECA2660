%A = load('data/data1.txt');
A = load('data/S3_128_100_diff_dataA.txt');
%B = load('data/data2.txt');
B = load('data/S3_128_100_diff_dataB.txt');
NTime = 64;
N = 128;
h = 1/N;
dt = h;
sigma = 1/32;% 2*h;
nu = 0;
c = 1;
Q = 1;
Ttime = NTime*dt;
gauss = @(x,t) Q./sqrt(pi*(sigma*sigma))* exp( -(mod(x,1.0)-c*t).^2 ./ (sigma*sigma)) ; 
%gauss = @(x,t) exp(-x.^2);
Npts = linspace(0,N-1,N-1);
Y = gauss((Npts-N/2)*h,NTime*dt);
plot(Npts,Y);
close all;
figure;
plot(A(:,1),A(:,2),A(:,1),A(:,3),'*');%,Npts,Y,'+');
%plot(Npts,Y);
%%
figure;
subplot(2,2,1);
loglog(B(:,1),B(:,2));
title('Qnh/Q');
subplot(2,2,2);
loglog(B(:,1),B(:,3)/B(1,3));
title('Enh/Enh(0)');
subplot(2,2,3);
loglog(B(:,1),B(:,4));
title('Rnh/sqrt(Enh(0)');
%%
   Rnh = 0;
   X   = A(:,1);
   U   = A(:,2);
   Uex = A(:,3);
   for j = 1 : length(X);
     Rnh = (U(j)-Uex(j))^2 + Rnh;
   end
   Rnh = sqrt(0.007812*Rnh)

%%
% kj = linspace(1,100,100);
% t = 0 ;
% 
% x = linspace(1,100,100);
% Uex = zeros(length(x),1);
% for n = 1 : length(x);
%     for j = 1 : length(kj)
%         Uex(n) = Uex(n) + exp( - kj(j).^2 * ( sigma.^2 )/4 ) * exp(i*kj(j)*(x(n)-c*t)); 
%     end    
% end
% plot(x,Uex*0.1);