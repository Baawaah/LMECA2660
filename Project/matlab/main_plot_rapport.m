%% Load

plot_all = 1;
plot_R   = 0;
plot_diag= 0;
plot_pres= 1;
%%
load('colormapsavefile.mat')
%t = [0 25 49 74 99 100];
t = [100];
%t =  [0 24];
%t = 99;
for k = 1 : length(t); 
    file_omega = sprintf('../data_server_3/CFD_omega_%d.txt',t(k));
    file_psi   = sprintf('../data_server_3/CFD_psi_%d.txt',t(k));
    file_u     = sprintf('../data_server_3/CFD_u_%d.txt',t(k));
    file_v     = sprintf('../data_server_3/CFD_v_%d.txt',t(k));
    file_R     = sprintf('../data_server_3/CFD_R_%d.txt',t(k));
    file_P     = sprintf('../data_server_3/CFD_P_%d.txt',0);
    file_R_pres= sprintf('../data_server_3/CFD_R_pres_%d.txt',0);
    omega (:,:,k) = load(file_omega);
    psi   (:,:,k) = load(file_psi);
    u     (:,:,k) = load(file_u);
    v     (:,:,k) = load(file_v);
    R     (:,:,k) = load(file_R);
    if plot_pres == 1
    P     (:,:,1) = load(file_P);
    R_pres(:,:,1) = load(file_R_pres);
    end
end
Diag = load('../data_server_3/CFD_DIAG.txt');
%% 
k = 1;
[SY SX] = size(omega(:,:,1));

% for l=1:SY
%     for m=1:SX
%         velocity(l,m,1) = sqrt(u(l,m,1)*u(l,m,1)+v(l,m,1)*v(l,m,1))/(1.5e-4);
%     end
% end
velocity = sqrt(u.*u+v.*v)/1.5e-4;
B = (u.*u+v.*v)/2+P(1:66,1:1000,:);

%%
if plot_all == 1 
for k = 1 : length(t)  
    
figure;
subplot(3,1,1);
im1 = imagesc(velocity(:,:,1));
axis equal; axis xy;
axis([0,SX,0,SY]);
caxis([min(min(velocity(:,:,k))) max(max(velocity(:,:,k)))]);
colormap(myColormap);
colorbar;
title('Dimensionless Velocity')
% s=sprintf('velocity_field_%d',k);
% print(s,'-depsc','-tiff')   

%figure;
subplot(3,1,2);
im2 = imagesc(psi(:,:,1)/1.5e-4);
axis equal; axis xy;
axis([0,SX,0,SY]);
caxis([min(min(psi(:,:,k)/1.5e-4)) max(max(psi(:,:,k)/1.5e-4))]);
colormap(myColormap);
colorbar;
title('Streamlines')
% s=sprintf('streamlines_%d',k);
% print(s,'-depsc','-tiff')

%figure;
subplot(3,1,3);
im3 = imagesc(omega(:,:,1)/1.5e-4);
axis equal; axis xy;
axis([0,SX,0,SY]);
caxis([min(min(omega(:,:,k)/1.5e-4)) max(max(omega(:,:,k)/1.5e-4))]);
colormap(myColormap);
colorbar;
title('Dimensionless Vorticity')
% s=sprintf('vorticity_%d',k);
% print(s,'-depsc','-tiff')
s=sprintf('plot_backward_tau2_%d',k);
print(s,'-depsc','-tiff')


figure;
subplot(2,1,1);
im4 = imagesc(P(:,:,1)/(1.5e-4*1.5e-4));
axis equal; axis xy;
axis([0,SX,0,SY]);
caxis([min(min(P(:,:,k)/(1.5e-4*1.5e-4))) max(max(P(:,:,k)/(1.5e-4*1.5e-4)))]);
colormap(myColormap);
colorbar;
title('Dimensionless Pressure')
% s=sprintf('pressure_%d',k);
% print(s,'-depsc','-tiff')

%figure;
subplot(2,1,2);
im5 = imagesc(B(:,:,1)/(1.5e-4*1.5e-4));
axis equal; axis xy;
axis([0,SX,0,SY]);
caxis([min(min(B(:,:,k)/(1.5e-4*1.5e-4))) max(max(B(:,:,k)/(1.5e-4*1.5e-4)))]);
colormap(myColormap);
colorbar;
title('Total Pressure')
% s=sprintf('pressure_tot_%d',k);
% print(s,'-depsc','-tiff')
s=sprintf('plot_backward_pressures_tau2_%d',k);
print(s,'-depsc','-tiff')

end    
end
%%
if plot_diag == 1
sample = [1 : 10 :length(Diag(:,1))];
figure;
RE_HL  = 10 * ones(length(sample));
RE_HOL = 5  * ones(length(sample));
plot(Diag(sample,1),Diag(sample,2),Diag(sample,1),Diag(sample,3),Diag(sample,1),RE_HL ,'--',Diag(sample,1),RE_HOL,'--');
%plot(Diag(:,1),Diag(:,3),Diag(:,1),RE_HOL,'--');
legend('Re_h','Re_h_\omega','Re_h Limit','Re_h_\omega Limit');
%legend('Re_h_\omega','Re_h_\omega Limit');
end    
  