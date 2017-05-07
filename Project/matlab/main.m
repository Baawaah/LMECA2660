close all;
t = 199;
file_omega = sprintf('../data/CFD_omega_%d.txt',t);
file_psi   = sprintf('../data/CFD_psi_%d.txt',t);
file_u     = sprintf('../data/CFD_u_%d.txt',t);
file_v     = sprintf('../data/CFD_v_%d.txt',t);
omega = load(file_omega);
psi   = load(file_psi);
u     = load(file_u);
v     = load(file_v);

%psi   = flipud(psi);
%omega = flipud(omega);
%u     = flipud(u);
%v     = flipud(v);

[Ni,Nj] = size(omega);
NY = linspace(1,Nj,Nj);
NX = linspace(1,Ni,Ni);
figure;

 hold on;
 %surf(NY,NX,omega,'edgecolor','none');
 %view(2)
% contour(NY,NX,omega,20);
 quiver(NY,NX,u,v);
 %contour(NY,NX,psi,20);
 %contour(NY,NX,omega,30);
 %caxis([0.002, 2])
 %colorbar
 hold off;