%% Load

plot_all = 1;
plot_R   = 0;
plot_diag= 1;
%%
t = [0 25 37 50 62 75 87 99];
%t = [0 25 50 75 99];
%t = 5;
for k = 1 : length(t); 
    file_omega = sprintf('../data/CFD_omega_%d.txt',t(k));
    file_psi   = sprintf('../data/CFD_psi_%d.txt',t(k));
    file_u     = sprintf('../data/CFD_u_%d.txt',t(k));
    file_v     = sprintf('../data/CFD_v_%d.txt',t(k));
    file_R     = sprintf('../data/CFD_R_%d.txt',t(k));
    omega(:,:,k) = load(file_omega);
    psi  (:,:,k) = load(file_psi);
    u    (:,:,k) = load(file_u);
    v    (:,:,k) = load(file_v);
    R    (:,:,k) = load(file_R);
end
Diag = load('../data/CFD_DIAG.txt');
%% 
k = 2;
[SY SX] = size(omega(:,:,1))
%%
if plot_all == 1 
for k = 1 : length(t)    
figure;
subplot(4,1,1);
im1 = imagesc(omega(:,:,k));
axis equal; axis xy; 
axis([0,SX,0,SY]);
caxis([-0.005 0.005]);
colormap jet;
colorbar;
title('Omega - Vorticity');
subplot(4,1,2);
im2 = contour(psi(:,:,k));
axis equal; axis xy;
axis([0,SX,0,SY]);
%caxis([-0.005 0.005]);
colormap jet;
colorbar;
title('Psi - Streamline');
subplot(4,1,3);
im3 = imagesc(u(:,:,k));
axis equal; axis xy;
axis([0,SX,0,SY]);
%caxis([-0.005 0.005]);
colormap jet;
colorbar;
title('u - Velocity')
subplot(4,1,4);
im4 = imagesc(v(:,:,k));
axis equal; axis xy;
axis([0,SX,0,SY]);
caxis([-0.005 0.005]);
colormap jet;
colorbar;
title('v - Velocity')
if plot_R == 1
figure;
im5 = imagesc(R(:,:,k));
axis equal; axis xy;
axis([0,SX,0,SY]);
caxis([-0.005 0.005]);
colormap jet;
colorbar;
title('R - Residu')
end    
end
end
%%
if plot_diag == 1
figure;
RE_HL  = 10 * ones(length(Diag(:,1)));
RE_HOL = 5  * ones(length(Diag(:,1)));
plot(Diag(:,1),Diag(:,2),Diag(:,1),Diag(:,3),Diag(:,1),RE_HL,'--',Diag(:,1),RE_HOL,'--');
legend('Re_h','Re_h_\omega','Re_h Limit','Re_h_\omega Limit');
end    
