%% Load

plot_all = 1;
plot_R   = 0;
plot_diag= 1;
plot_pres= 0;
%%
load('colormapsavefile.mat')
%t = [0 25 37 50 62 75 87 99];
t = [0 50 100];
%t =  [0 24];
%t = 99;
for k = 1 : length(t); 
    file_omega = sprintf('../data/CFD_omega_%d.txt',t(k));
    file_psi   = sprintf('../data/CFD_psi_%d.txt',t(k));
    file_u     = sprintf('../data/CFD_u_%d.txt',t(k));
    file_v     = sprintf('../data/CFD_v_%d.txt',t(k));
    file_R     = sprintf('../data/CFD_R_%d.txt',t(k));
    file_P     = sprintf('../data/CFD_P_%d.txt',t(k));
    file_R_pres= sprintf('../data/CFD_R_pres_%d.txt',t(k));
    omega (:,:,k) = load(file_omega);
    psi   (:,:,k) = load(file_psi);
    u     (:,:,k) = load(file_u);
    v     (:,:,k) = load(file_v);
    R     (:,:,k) = load(file_R);
    if plot_pres == 1
    P     (:,:,k) = load(file_P);
    R_pres(:,:,k) = load(file_R_pres);
    end
end
Diag = load('../data/CFD_DIAG.txt');
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
sample = [1 : 10 :length(Diag(:,1))];
figure;
RE_HL  = 10 * ones(length(sample));
RE_HOL = 5  * ones(length(sample));
plot(Diag(sample,1),Diag(sample,2),Diag(sample,1),Diag(sample,3),Diag(sample,1),RE_HL ,'--',Diag(sample,1),RE_HOL,'--');
%plot(Diag(:,1),Diag(:,3),Diag(:,1),RE_HOL,'--');
legend('Re_h','Re_h_\omega','Re_h Limit','Re_h_\omega Limit');
%legend('Re_h_\omega','Re_h_\omega Limit');
end    
%% 
if plot_pres == 1
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