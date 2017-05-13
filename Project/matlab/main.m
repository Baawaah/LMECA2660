%close all;
t = 19999;
file_omega = sprintf('../data/CFD_omega_%d.txt',t);
file_psi   = sprintf('../data/CFD_psi_%d.txt',t);
file_u     = sprintf('../data/CFD_u_%d.txt',t);
file_v     = sprintf('../data/CFD_v_%d.txt',t);
file_R     = sprintf('../data/CFD_R_%d.txt',t);
omega = load(file_omega);
psi   = load(file_psi);
u     = load(file_u);
v     = load(file_v);
R     = load(file_R);

%psi   = flipud(psi);
%omega = flipud(omega);
%u     = flipud(u);
%v     = flipud(v);

[Ni,Nj] = size(omega);
NY = linspace(1,Nj,Nj); 
NX = linspace(1,Ni,Ni);
%%
figure;
subplot(2,2,1);
quiver(NY,NX,u,v);
subplot(2,2,2);
contour(NY,NX,R,50);
colorbar
subplot(2,2,3);
contour(NY,NX,psi,50);
subplot(2,2,4);
contour(NY,NX,omega,50);
colorbar
%%
figure;
im1 = imagesc(omega);
axis equal; axis xy;
caxis([-0.005 0.005]);
colormap jet;
colorbar;
%%
 hold on;
 %surf(NY,NX,omega,'edgecolor','none');
 %view(2)
 %contour(NY,NX,omega,20);
 %quiver(NY,NX,u,v);
 %contour(NY,NX,u,100);
 %contour(NY,NX,psi,50);
 %contour(NY,NX,R,100);
 %contour(NY,NX,omega,100);
 %caxis([0.002, 2])
 %colorbar
 hold off;

 figure;
 surf(u,'edgecolor','none');
 view(2);
 colormap jet;
 colorbar;
 %%