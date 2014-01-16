load wannier_FeSe_4d_matrix_v2
load spatial_ldos2_FeSe_N_15_U_0955_Vimp_50_E_-0.0084_ita_0.001_z_0.mat
zGridPoint = 0;
E = .0084;
localLdos = ldos;
a = 7.23;
z = zGrid(41 + zGridPoint);
[X, Y] = meshgrid(xGridRange,yGridRange);
h = figure;
pcolor(X,Y,localLdos');
axis('square')
colorbar;
titleName = ['E = ', num2str(E), ', z = ', num2str(z), ' Bohr'];
title(titleName);
xlabel('x (Bohr)');
ylabel('y (Bohr)');

