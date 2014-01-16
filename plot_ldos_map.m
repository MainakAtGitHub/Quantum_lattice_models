function p=plot_ldos_map(ldosfile)
load wannier_FeSe_4d_matrix_v2
load(ldosfile,'-mat')
zGridPoint = 0;
E = .0084;
localLdos = loacalLdos;
a = 7.23;
z = zGrid(41 + zGridPoint);
[X, Y] = meshgrid(xGridRange,yGridRange);
% Create figure
figure1 = figure;
% Create axes
%axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1]);
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'CameraViewAngle',3.88994451795861);
xlim(axes1,[min(X(:)) max(X(:))]);
ylim(axes1,[min(Y(:)) max(Y(:))]);
grid(axes1,'on');
hold(axes1,'all');
surf(X,Y,localLdos','LineStyle','none','FaceColor','interp');
% view from 1 Fe zone!
view(axes1,[45 90]);
%view([0 90])
%pcolor(X,Y,localLdos');
%axis('square')
colorbar;
titleName = ['E = ', num2str(E), ', z = ', num2str(z), ' Bohr'];
title(titleName);
xlabel('x (Bohr)');
ylabel('y (Bohr)');

