%% Load

plot_all = 0;
plot_R   = 0;
plot_diag= 1;
plot_pres= 0;
plot_fin = 1;
%%
load('colormapsavefile.mat')
%t = [0 25 37 50 62 75 87 99];
%t = [0 50 99 149 200];
%t =  [0 24];
%t = 99;
%t = 100;
%t = [ 5 24 34 54 74 84 100]
%t= [66 70 73 76 100]
t = 100
for k = 1 : length(t); 
%     file_omega = sprintf('../data_server_11_t100_R150_h0015/CFD_omega_%d.txt',t(k));
%     file_psi   = sprintf('../data_server_11_t100_R150_h0015/CFD_psi_%d.txt',t(k));
%     file_u     = sprintf('../data_server_11_t100_R150_h0015/CFD_u_%d.txt',t(k));
%     file_v     = sprintf('../data_server_11_t100_R150_h0015/CFD_v_%d.txt',t(k));
%     file_R     = sprintf('../data_server_11_t100_R150_h0015/CFD_R_%d.txt',t(k));
%     file_P     = sprintf('../data_server_11_t100_R150_h0015/CFD_P_%d.txt',t(k));
%     file_R_pres= sprintf('../data_server_11_t100_R150_h0015/CFD_R_pres_%d.txt',t(k));

    file_omega = sprintf('../data_pc_2_t100_R150_h0015_f/CFD_omega_%d.txt',t(k));
    file_psi   = sprintf('../data_pc_2_t100_R150_h0015_f/CFD_psi_%d.txt',t(k));
    file_u     = sprintf('../data_pc_2_t100_R150_h0015_f/CFD_u_%d.txt',t(k));
    file_v     = sprintf('../data_pc_2_t100_R150_h0015_f/CFD_v_%d.txt',t(k));
    file_R     = sprintf('../data_pc_2_t100_R150_h0015_f/CFD_R_%d.txt',t(k));
    file_P     = sprintf('../data_pc_2_t100_R150_h0015_f/CFD_P_%d.txt',t(k));
    file_R_pres= sprintf('../data_pc_2_t100_R150_h0015_f/CFD_R_pres_%d.txt',t(k));


%     file_omega = sprintf('../data/CFD_omega_%d.txt',t(k));
%     file_psi   = sprintf('../data/CFD_psi_%d.txt',t(k));
%     file_u     = sprintf('../data/CFD_u_%d.txt',t(k));
%     file_v     = sprintf('../data/CFD_v_%d.txt',t(k));
%     file_R     = sprintf('../data/CFD_R_%d.txt',t(k));
%     file_P     = sprintf('../data/CFD_P_%d.txt',t(k));
%     file_R_pres= sprintf('../data/CFD_R_pres_%d.txt',t(k));

    omega (:,:,k) = load(file_omega);
    psi   (:,:,k) = load(file_psi);
    u     (:,:,k) = load(file_u);
    v     (:,:,k) = load(file_v);
    R     (:,:,k) = load(file_R);
    P     (:,:,k) = load(file_P);
    R_pres(:,:,k) = load(file_R_pres);

end
% Diag = load('../data_server_11_t100_R150_h0015/CFD_DIAG.txt');
% Diag = load('../data/CFD_DIAG.txt');
 Diag = load('../data_pc_2_t100_R150_h0015_f/CFD_DIAG.txt');
%% 
k = 1;
[SY SX] = size(omega(:,:,1));
%%
if plot_all == 1 
for k = 1 : length(t)    
figure;
subplot(4,1,1);
im1 = imagesc(omega(:,:,k));
axis equal; axis xy; 
axis([0,SX,0,SY]);
caxis([-0.003 0.003]);
colormap(myColormap);
colorbar;
title('Omega - Vorticity');
subplot(4,1,2);
im2 = contour(psi(:,:,k));
axis equal; axis xy;
axis([0,SX,0,SY]);
%caxis([-0.005 0.005]);
colormap(jet);
colorbar;
title('Psi - Streamline');
subplot(4,1,3);
im3 = imagesc(u(:,:,k));
axis equal; axis xy;
axis([0,SX,0,SY]);
caxis([-0.0003 0.0003]);
colormap(myColormap);
colorbar;
title('u - Velocity')
subplot(4,1,4);
im4 = imagesc(v(:,:,k));
axis equal; axis xy;
axis([0,SX,0,SY]);
caxis([-0.0002 0.0002]);
colormap(myColormap);
colorbar;
title('v - Velocity')
if plot_R == 1
figure;
im5 = imagesc(R(:,:,k));
axis equal; axis xy;
axis([0,SX,0,SY]);
%caxis([-0.0001 0.0001]);
colormap(myColormap);
colorbar;
title('R - Residu')
end    
end
end
%%
if plot_diag == 1
sample = [1 : 600 :length(Diag(:,1))];
figure;
RE_HL  = 10 * ones(length(sample));
RE_HOL = 5  * ones(length(sample));
%hold on;
plot(Diag(sample,1),Diag(sample,2),Diag(sample,1),Diag(sample,3),Diag(sample,1),RE_HL,'--',Diag(sample,1),RE_HOL,'--');
%plot(Diag(sample,1),RE_HL,Diag(sample,1),RE_HOL);
%hold off;
%plot(Diag(:,1),Diag(:,3),Diag(:,1),RE_HOL,'--');
legend('Re_h','Re_h_\omega','Re_h Limit','Re_h_\omega Limit');

title('Diagnostic Re_h and Re_h_\omega')
xlabel('\tau - Adimensional time')
axis([0,1,0,11]);
%legend('Re_h_\omega','Re_h_\omega Limit');
figure;
plot(Diag(sample,1),Diag(sample,4));
end    
%% 
if plot_pres == 1
k = 1
figure;
imP1 = imagesc(P(:,:,k));
axis equal; axis xy;
axis([0,SX,0,SY]);
%caxis([-0.0001 0.0001]);
colormap(myColormap);
colorbar;
title('P - Pressure')
figure;
imP2 = imagesc(R_pres(:,:,k));
axis equal; axis xy;
axis([0,SX,0,SY]);
%caxis([-0.0001 0.0001]);
colormap(myColormap);
colorbar;
title('R - Residu')    
end    
%%
% figure;
% imP2 = imagesc(omega(:,:,2));
% axis equal; axis xy;
% axis([0,SX,0,SY]);
% caxis([-0.001 0.001])
% colormap(myColormap);
%%
if plot_fin == 1
k = 1;

figure;
UV = sqrt(u(:,:,k).*u(:,:,k) + v(:,:,k).*v(:,:,k))/(1.5e-4);
im1 = imagesc(linspace(0,15,SX),linspace(0,1,SY),UV(:,:,k));
xlabel('x/H');
ylabel('y/H');
axis equal; axis xy;
axis([0,15,0,1]);
colormap(myColormap);
c1 = colorbar('eastoutside'); 
set(c1, 'Position',[0.92    0.4845    0.02    0.07]);
title('Dimensionless Velocity Field');
%%
figure;
omegaQ = omega/(1.5e-4);
im4 = imagesc(linspace(0,15,SX),linspace(0,1,SY),omegaQ(:,:,k));
xlabel('x/H');
ylabel('y/H');
axis equal; axis xy;
axis([4,6,0,1]);
caxis([-30 300])
colormap(myColormap);
c4 = colorbar('eastoutside'); 
set(c4, 'Position',[0.92    0.4845    0.02    0.07]);
title('Dimensionless Vorticity Field');
%%
figure;
psiQ = psi/(1.5e-4);
im2 = contour(linspace(0,15,SX),linspace(0,1,SY),psiQ(:,:,k),25);
xlabel('x/H');
ylabel('y/H');
axis equal; axis xy;
axis([4,6,0,1]);
colormap(myColormap);
c2 = colorbar('eastoutside'); 
set(c2, 'Position',[0.92    0.4845    0.02    0.07]);
title('Dimensionless Streamline');
%%
figure;
k = 1;
PHQ = (P-P(SY/2-1,SX+1,k))/((1.5e-4)^2);
im2 = imagesc(linspace(0,15,SX+1),linspace(0,1,SY+1),PHQ(:,:,k));
xlabel('x/H');
ylabel('y/H');
axis equal; axis xy;
axis([0,15,0,1]);
colormap(myColormap);
c3 = colorbar('eastoutside'); 
set(c3, 'Position',[0.92    0.4845    0.02    0.07]);
title('Dimensionless Pressure Field');

% figure;
% UV = sqrt(u(:,:,k).*u(:,:,k) + v(:,:,k).*v(:,:,k))/(1.5e-4);
% x = linspace(0,SX,SX);
% y = linspace(0,SY,SY);
% imagesc(x,y,UV);
% 
% axis equal; axis xy;
% axis([0,SX,0,SY]);
% xticklabels = 0:20:100;
% caxis([0.0 1]);
% colormap(myColormap);
% colorbar;
% title('Velocity Field');
% 
% 
% L=15;
% H=1;
% Ls=5*H;
% Hs=1/2;
% for k = int64([H (Ls-H) (Ls-H/4) Ls (Ls+H/4) (Ls+H) (L-H)]*(1/0.015))  
% figure
% plot(abs(u(:,k,1)/1.5e-4),linspace(0,1,SY))
% title('Dimensionless Horizontal Velocity Profile')
% ylabel('y/H')
% xlabel('uH/Q_0')
% % s=sprintf('adim_velocity_field_%d',k);
% % print(s,'-depsc','-tiff')   
% % figure
% % plot(abs(omega(:,k,1)/1.5e-4),linspace(0,1,SY))
% % title('Dimensionless Horizontal Velocity Profile')
% % ylabel('y/H')
% % xlabel('\omega H^2/Q_0')
% % % s=sprintf('adim_vorticity_field_%d',k);
% % % print(s,'-depsc','-tiff')  
% 
% figure;
% %subplot(2,1,1);
% im4 = imagesc(P(:,:,1)/(1.5e-4*1.5e-4));
% axis equal; axis xy;
% axis([0,SX,0,SY]);
% caxis([min(min(P(:,:,k)/(1.5e-4*1.5e-4))) max(max(P(:,:,k)/(1.5e-4*1.5e-4)))]);
% colormap(myColormap);
% colorbar;
% title('Dimensionless Pressure')
% % s=sprintf('pressure_%d',k);
% % print(s,'-depsc','-tiff')
% figure;
% %subplot(2,1,2);
% im5 = imagesc(B(:,:,1)/(1.5e-4*1.5e-4));
% axis equal; axis xy;
% axis([0,SX,0,SY]);
% caxis([min(min(B(:,:,k)/(1.5e-4*1.5e-4))) max(max(B(:,:,k)/(1.5e-4*1.5e-4)))]);
% colormap(myColormap);
% colorbar;
% title('Total Pressure')
% % s=sprintf('pressure_tot_%d',k);
% % print(s,'-depsc','-tiff')
% s=sprintf('plot_backward_pressures_tau2_%d',k);
% print(s,'-depsc','-tiff')
% end



end

%%