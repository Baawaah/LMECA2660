omega = load('../data/CFD_omega_40.txt');
psi   = load('../data/CFD_psi_40.txt');

psi   = flipud(omega);
omega = flipud(omega);

[Ni,Nj] = size(omega);
NY = linspace(1,Nj,Nj);
NX = linspace(1,Ni,Ni);
hold on;
%mesh(NY,NX,omega),'edgecolor','none');
%view(2)
contour(NY,NX,omega,20);
hold off;