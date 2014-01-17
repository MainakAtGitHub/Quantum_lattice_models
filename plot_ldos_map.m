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
% range=3
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'CameraViewAngle',3.88994451795861);
%range=5
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'CameraViewAngle',1.94553273895407);
xlim(axes1,[min(X(:)) max(X(:))]);
ylim(axes1,[min(Y(:)) max(Y(:))]);
grid(axes1,'on');
hold(axes1,'all');
surf(X,Y,localLdos','LineStyle','none','FaceColor','interp');
% view from 1 Fe zone!
view(axes1,[45 90]);
 bluemap(figure1)
%view([0 90])
%pcolor(X,Y,localLdos');
%axis('square')
colorbar;
titleName = ['E = ', num2str(E), ', z = ', num2str(z), ' Bohr'];
title(titleName);
xlabel('x (Bohr)');
ylabel('y (Bohr)');
thickness=0.4;
ratio=0.8;
lgth=12;
offset=[-10,10];
x=[0,lgth]+offset(1);
y=[0,lgth]+offset(2);
z=[1 1]*max(localLdos(:));
h1=arrow3d(x,y,z,ratio,thickness);
text(1.1*x(2),1.1*y(2),z(2),'$x$','FontSize',24,'Interpreter','latex')
set(h1,'facecolor',[1 0 0])
x=[0,-lgth]+offset(1);
y=[0,lgth]+offset(2);
z=[1 1]*max(localLdos(:));
h1=arrow3d(x,y,z,ratio,thickness);
text(1.1*x(2),1.1*y(2),z(2),'$y$','FontSize',24,'Interpreter','latex')
% label the coordinate axes
set(h1,'facecolor',[1 0 0])
if isunix
    % create pdf of figure
    [~,filename,extension]=fileparts(ldosfile);
    print('-djpeg', ['/tmp/',filename,extension,'.jpg']);
%print_pdf(['/tmp/',filename,extension,'cut',num2str(N),'.pdf'])
end;
